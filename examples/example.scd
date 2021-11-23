(
s.options.memSize = 8192 * 64;
s.boot
)


// SOUND
// SEED
// AUTOMATA


(
// init
~ssa = SoundSeedAutomata.new(
    s, SoundSeedGarden.new(s));
)

(
// SoundSeedGarden eats a lot of CPU
// can also use a bus number,
// skipping the dynamics and ambience
// ~ssa = SoundSeedAutomata.new(s, 0);
)


(
//plant a seed
var tag;
Routine {
    tag = ~ssa.plant;
    0.5.wait;
    ~ssa.sprout(tag);
}.play
)


//
// ~ssa.garden.send_bus.scope
// ~ssa.garden.sum_bus.scope
// ~ssa.garden.exp_bus.scope
// ~ssa.garden.max_bus.scope
//
// ~ssa.garden.mixdown_synth.trace
//
// ~ssa.garden.ambience_synth.trace
//
// ~ssa.seeds[1].synth.trace
//
// ~ssa.seeds[1].soundbuf.plot
//
//
// ~ssa.seeds[6].free
//
// s.queryAllNodes
//
// OSCFunc.trace


s.record
s.stopRecording
s.quit