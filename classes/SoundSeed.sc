SoundSeedCrossover {
    *ar {arg sig, freq;
        var lp = LPF.ar(LPF.ar(sig,freq),freq), hp = HPF.ar(HPF.ar(sig,freq),freq);
        // var lp = LPF.ar(sig,freq,mul:-1), hp = HPF.ar(sig,freq);
        ^[lp, hp]
    }
}

SoundSeedCrossoverTree {
    // tree of crossover filters
    *ar {
        arg input,
        lo = 100,
        hi = 8000,
        log_n = 3;

        var signals = [input];
        var freqs = [(lo.log + hi.log / 2).exp];
        var base = hi/lo ** ( 1/ (log_n-1*log_n));
        log_n.do{
            arg level;
            var spread = log_n - 1 - level;
            spread;
            signals = signals.collect{arg sig, idx;
                SoundSeedCrossover.ar(sig, freqs[idx])}.flatten;
            freqs = freqs.collect{arg freq;
                [base**(0-spread)*freq, base**spread*freq]}.flatten
        };
        ^signals
    }
}

SoundSeedPan {
    //pseudo-binaural pan
    *ar {
        arg signal, pan=0;
        var path_dist = 0.3;
        var max_dist = 10;
        var head_radius = 0.1;
        var sound_speed = 300;
        var min_time = SampleDur.ir * 4.0;
        var disps = [0-head_radius, head_radius]-pan;
        var dists = (path_dist*path_dist + (disps*disps)).sqrt;
        var gains = 1/(1+dists); // inverse square instead?
        var signals = [signal, signal]*0.71;
        var max_time = max_dist/sound_speed+min_time;
        signals = BLowShelf.ar(signals, 1200, db:gains.log10*20);
        ^DelayC.ar(signals, max_time,
            (dists/sound_speed+min_time).min(max_time),
            mul:gains);
    }
}

SoundSeed {
    // TODO: get rid of databuf, pass length back into synth?
    var <server;
    var <inbus;

    const numChannelsIn = 1;
    var max_frames;

    var <cond, <rec_length, <recorder, <synth;
    var <soundbuf, <databuf, <featbuf, <statsbuf;

    var <state_bits = 5;
    var <symbol_bits = 4;
    var <move_bits = 1;
    var <step = 0;
    var <head = 0;
    var <state = 0;
    var <movements = #[-1, 1];
    var <program_length;
    var <program_width;
    var <program_bits;
    var <program_bytes;
    var <tape;
    var <visits;
    var <program;

    var <intervals;

    var mel_bands = 64;
    var n_moments = 4; // max 4 (7)
    var n_derivs = 2; //  max 2

    *new { arg server, inbus;
        ^super.newCopyArgs(server, inbus).init;
    }

    init {
        max_frames = server.sampleRate * 10;
        soundbuf = Buffer.alloc(server, max_frames, numChannelsIn);
        databuf = Buffer.alloc(server, 1);

        featbuf = Buffer.new;
        statsbuf = Buffer.new;

        cond = CondVar.new;

        intervals = [1, 3/2, 5/4, 2, 4/3, 5/3, 9/8, 15/8];

        program_length = 1<<(state_bits+symbol_bits);
        program_width = state_bits+symbol_bits+move_bits;
        program_bits = (program_length*program_width);
        program_bytes = (program_bits/8).ceil.div(1);

        tape = Dictionary.new;
        visits = Dictionary.new;

        recorder = Synth(\sound_seed_recorder, [
            inbus:inbus, soundbuf:soundbuf, databuf:databuf
        ]);

        OSCFunc({ arg msg, time;
            var path, node_id, trig_id, val;
            # path, node_id, trig_id, val = msg.postln;
            rec_length = val;
            cond.signalOne;
        },'/tr', server.addr, argTemplate:[recorder.nodeID.postln]).oneShot;
    }

    // germinate {}

    read_offset { arg byte, width, offset;
        var err = if(byte>=program_bytes) {
            "program out of bounds at %".format(byte).postln};
        var b0 = program[byte];
        var b1 = program[byte+1]; //program storage should have an extra byte?
        var mask = (0xff >> (8 - width));
        ^(
            ((b0 & 0xff) << (offset + width - 8))
            | ((b1 & 0xff) >> (16 - offset - width))
        ) & mask
    }

    get_instruction { arg state, symbol;
        var idx = (state<<symbol_bits) | symbol;
        var state_bit = idx*program_width;
        var state_byte = state_bit.div(8);
        var state_offset = state_bit - (state_byte*8);
        var symbol_bit = state_bit + state_bits;
        var symbol_byte = symbol_bit.div(8);
        var symbol_offset = symbol_bit - (symbol_byte*8);
        var move_bit = symbol_bit + symbol_bits;
        var move_byte = move_bit.div(8);
        var move_offset = move_bit - (move_byte*8);
        var new_state  = this.read_offset(state_byte,  state_bits,  state_offset);
        var new_symbol = this.read_offset(symbol_byte, symbol_bits, symbol_offset);
        var move_idx   = this.read_offset(move_byte,   move_bits,   move_offset);
        ^[new_state, new_symbol, move_idx]
    }

    get_program {
        var quantize = {
            arg val, levels, floor = -13, range = 10;
            ((val - floor).clip(0, range-1e-7)/range*levels).div(1)
        };

        program = Int8Array.newClear(program_bits.div(8)+4);

        //start of audio -> program
        /*soundbuf.loadToFloatArray(0, program_bits.div(32)+1, {
        arg arr;
        arr.do{ arg val, idx;
        var bin = val.as32Bits;//.postln;
        var pidx = idx*4;
        program[pidx] = (bin >> 24) & 0xff;
        program[pidx+1] = (bin >> 16) & 0xff;
        program[pidx+2] = (bin >> 8) & 0xff;
        program[pidx+3] = bin & 0xff;
        };
        });*/

        // multichannel buffers come back channel-first (frame locality)

        //mfcc stats -> program

        //quantize to program_width so 1 stat per instruction
        //this would work for a "sparse" Int32 program format:
        // "% stats vs % instructions".format(n_moments*(n_derivs+1)*mel_bands, program_length).postln;
        /*statsbuf.loadToFloatArray(0, program_length, {
        arg arr;
        var pidx = 0;
        var levels = 2**program_width;
        (n_derivs+1).do{
        arg deriv;
        7.do{
        arg moment;
        if (moment < n_moments){mel_bands.do{
        arg coef;
        var idx = coef + (mel_bands*(moment + 7*deriv));
        var stat = arr[idx];
        program[pidx] = quantize.(stat, levels, floor:-20, range:60).postln;
        pidx = pidx+1;
        }}
        }
        }
        });*/
        "% stats vs % bytes".format(
            n_moments*(n_derivs+1)*mel_bands, program_bits.div(8)).postln;
        // (statsbuf.numFrames*statsbuf.numChannels).postln;
        statsbuf.loadToFloatArray(0, program_length, {
            arg arr;
            var pidx = 0;
            (n_derivs+1).do{
                arg deriv;
                7.do{
                    arg moment;
                    if (moment < n_moments){mel_bands.do{
                        arg coef;
                        var idx = coef + (mel_bands*(moment + 7*deriv));
                        var stat = arr[idx];
                        if(pidx < program_bits.div(8)){
                            program[pidx] = quantize.(
                                stat.abs, 256, floor:0, range:60).bitXor(0xf0);
                        };
                        pidx = pidx+1;
                    }}
                }
            }
        });
    }

    sprout {
        synth = Synth(\sound_seed_default_synth, [
            soundbuf:soundbuf, databuf:databuf
        ]);
        Routine {
            recorder.set(\stop, 1);
            cond.wait;
            recorder.free;

            \stopped.postln;

            FluidBufMFCC.process(server, soundbuf, 0, rec_length, features:featbuf,
                numBands:mel_bands, numCoeffs:mel_bands,
                action:{cond.signalOne});
            cond.wait;

            FluidBufStats.process(server, featbuf, 0, featbuf.numFrames, stats:statsbuf,
                numDerivs:n_derivs,
                action:{cond.signalOne});
            cond.wait;

            /*FluidBufChroma.process(s, soundbuf, 0, rec_length, 0, -2, chromabuf,
                numChroma:n_chroma, ref:440, minFreq:40, maxFreq: 4096,
                windowSize:chroma_window, hopSize:chroma_window.div(2), fftSize:4096,
                // normalize:1,
                action:{~node_cond[node_id].signalOne});
            ~node_cond[node_id].wait;*/


            \features.postln;

            this.get_program.();

            \program.postln;

            loop {
                var vis, vis_key;
                var symbol = tape[head] ? 0;
                var new_state, new_symbol, move_idx, move;
                // "step: % \t head: % \t state: % \t symbol: %".format(step, head, state, symbol).postln;
                # new_state, new_symbol, move_idx = this.get_instruction(
                    state, symbol);

                // tricks to get more interesting machines:
                // move = if(symbol==0){0-head.sign}{movements[move_idx]}; //symbol zero move -> center
                // new_symbol = if(new_symbol==0){1}{new_symbol}; // 0 is never written
                move = movements[move_idx];

                if (new_symbol!=0) {tape[head] = new_symbol};
                head = head + move;

                vis_key = [head, state, symbol];
                vis = (visits[vis_key] ? 0) + 1;
                visits[vis_key] = vis;

                state = new_state;
                step = step + 1;
                synth.set(
                    \rate, intervals[symbol>>1] * (2**(symbol&1)) / 2,
                    \center, head/100,
                    \distance, head.abs/100 + (vis/10),
                );
                ((state>>1 + 1) * 0.03 / if((state&1)==0){1.5}{1}).wait;
            };
        }.play;
    }
}

SoundSeedGarden {
    var <server;
    var <log_n_bands;

    var <n_bands;
    var <reduce;
    var <cond;
    var <outbus;

    var <send_bus, <sum_bus, <exp_bus, <max_bus;
    var <mixdown_synth, <ambience_synth;

    *new { arg server, log_n_bands=3, outbus=0;
        ^super.newCopyArgs(server, log_n_bands, outbus).init;
    }

    init {
        n_bands = 1<<log_n_bands;
        send_bus = Bus.audio(server, 2); // raw, single-band, stereo, mixdown
        sum_bus  = Bus.audio(server, n_bands*2); // raw, multi-band, stereo, mixdown
        exp_bus  = Bus.audio(server, n_bands*2); // expanded, multi-band, stereo mixdown
        max_bus  = Bus.audio(server, n_bands); // envelope, multi-band, mono, reduce-max

        reduce = {
            arg input, distance, sharpness=3,
            rms_cutoff=100, attack=3e-3, release=1e-1;
            var delay = 8*SampleDur.ir + attack;
            var edist = 2**distance;
            var distanced = LPF.ar(HPF.ar(input,
                (20*edist).min(20000)), (20000/edist).max(40))/edist;
            var to_send = input/(distance+1);
            var bands_l = SoundSeedCrossoverTree.ar(distanced[0], log_n:log_n_bands);
            var bands_r = SoundSeedCrossoverTree.ar(distanced[1], log_n:log_n_bands);
            var bands = bands_l+bands_r;
            var envs = Amplitude.ar(
                BLowPass4.ar(bands*bands, rms_cutoff).sqrt, attack, release);
            var delayed = DelayN.ar(bands_l++bands_r, 0.1, delay);
            var exp_envs = (envs*sharpness).exp;
            var expanded = delayed * exp_envs;
            // distanced.poll(10, \distanced); input.poll(10, \input); envs.poll(10, \env);
            [
                Out.ar(send_bus, to_send), // singleband stereo
                Out.ar(sum_bus, delayed), // multiband stereo
                Out.ar(exp_bus, expanded), // multiband stereo
                Out.ar(max_bus, (envs-In.ar(max_bus, n_bands)).max(0)) //multiband mono
            ]
        };

        SynthDef(\sound_seed_mixdown, { arg outbus=0;
            var rms_cutoff = \rms_cutoff.ir(100);
            var max_thresh = \max_thresh.ir(0.2);
            var min_thresh = \min_thresh.ir(0.01);
            var ratio = \ratio.ir(0.125);
            var attack = \attack.ir(3e-3);
            var release = \release.ir(1e-1);

            var delay = 8*SampleDur.ir + attack;
            var bands_raw_l = In.ar(sum_bus, n_bands);
            var bands_raw_r = In.ar(sum_bus.index+n_bands, n_bands);
            var bands_raw = bands_raw_l + bands_raw_r;
            var bands_exp_l = In.ar(exp_bus, n_bands);
            var bands_exp_r = In.ar(exp_bus.index+n_bands, n_bands);
            var bands_exp = bands_exp_l + bands_exp_r;
            // filter after max here: is delay still right?
            // var env_max = BLowPass4.ar(In.ar(max_bus, n_bands), rms_cutoff);
            var env_max = In.ar(max_bus, n_bands);

            var env_raw = Amplitude.ar(
                BLowPass4.ar(bands_raw*bands_raw, rms_cutoff).sqrt.sanitize,
                attack, release);
            var env_exp = Amplitude.ar(
                BLowPass4.ar(bands_exp*bands_exp, rms_cutoff).sqrt.sanitize,
                attack, release);
            var thresh = env_max.min(max_thresh).max(min_thresh);
            var env = (env_raw - thresh).max(0) * ratio + env_raw.min(thresh); //mono
            var bands_exp_delayed_l = DelayN.ar(bands_exp_l, 0.1, delay);
            var bands_exp_delayed_r = DelayN.ar(bands_exp_r, 0.1, delay);
            var gain = env / (env_exp+min_thresh);
            var bands_out_l = bands_exp_delayed_l * gain;
            var bands_out_r = bands_exp_delayed_r * gain;

            var sig_out = [bands_out_l.sum, bands_out_r.sum]; // over bands
            // env_raw.poll(10, \raw); env_exp.poll(10, \exp);

            Out.ar(outbus, sig_out);
        }).add;

        // this seems fragile? how can I know it's the right done message?
        OSCFunc({ arg msg, time;
            ("This is the done message for the SynthDef.send:" + [time, msg]).postln;
            mixdown_synth = Synth.tail(server, \sound_seed_mixdown, [\outbus, outbus]);
        }, '/done').oneShot; // remove me automatically when done

        // TODO: default GVerb here, external option
        ambience_synth = {
            MiVerb.ar(In.ar(send_bus, 2), mul:0.05, drywet:1, time:0.1)
        }.play(addAction:'addToTail');
    }
}

SoundSeedAutomata {
    var <server;
    var <garden;
    var <default_inbus;

    const <numChannelsIn = 1;
    const <numChannelsOut = 2;
    var <seeds;
    var <counter = 0;
    var last_tag;

    *new { arg server, garden=nil, default_inbus=0;
        ^super.newCopyArgs(server, garden, default_inbus.postln).init;
    }

    init {
        seeds = Dictionary.new;

        SynthDef(\sound_seed_recorder, { arg inbus, soundbuf, databuf;
            var input = SoundIn.ar(inbus, numChannelsIn);
            var stop = T2A.ar(\stop.tr(0));
            var phase = Phasor.ar(
                1, BufRateScale.kr(soundbuf), start:0, end:BufFrames.kr(soundbuf));
            Demand.ar(stop, 0, Dbufwr(phase, databuf, 0));
            SendTrig.ar(stop, 0, phase);
            BufWr.ar(input, soundbuf, phase, loop: 0);
            // input.poll(10, \rec)
        }).add;

        SynthDef(\sound_seed_default_synth, { arg soundbuf, databuf;
            //  if there is a garden, ignore outbus and use garden's `reduce`
            var distance = \distance.kr(0, 0.1);
            var output = if (garden.class==Integer) {
                Out.ar(garden, _)
            }{
                // {arg sig; garden.reduce(sig, distance:distance)}
                garden.reduce.(_, distance:distance)
            };
            var frames = BufRd.ar(numChannelsIn, databuf, DC.ar(0), interpolation:1);
            var phase = Phasor.ar(0,
                (frames>0)*\rate.kr(0, 0.02)*BufRateScale.kr(soundbuf), 0, frames);
            var sound = BufRd.ar(numChannelsIn, soundbuf, phase);
            var panned = SoundSeedPan.ar(
                sound * \gain.kr(0.5, 0.3),
                \center.kr(0, 0.03)
            );
            // sound.poll(10, \sound); panned.poll(10, \panned);
            output.(panned)
        }).add;
    }

    plant { arg tag=nil, inbus=nil;
        inbus = inbus ? default_inbus;
        last_tag = tag = tag ?? {counter = counter + 1; counter};
        seeds[tag] = SoundSeed.new(server, inbus);
        "planted seed '%'".postf(tag);
        ^tag
    }

    sprout { arg tag;
        tag = tag ? last_tag;
        ^seeds[tag].sprout;
    }
}