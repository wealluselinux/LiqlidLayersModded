# Liquid Layers Modded

A from-scratch rebuild of the "Liquid Layers" style layered-liquid simulator, with a mod: drop objects into the tank and they float, sink or get flung around by the liquid.

Inspired by [Liquid Layers](https://grantkot.com/ll/) by Grant Kot. This is an independent reimplementation, not a copy of that site's code.

## Objects

| Object | Behaviour |
|--------|-----------|
| Ball | Light, bobs on the top layer and spins when the liquid moves |
| Duck | Light, weighted low so it rights itself |
| Crate | Medium density, settles partway down between layers |
| Boat | Keel-weighted hull with a sail, stays upright |
| Anchor | Dense, sinks to the bottom |

## Controls

- Pick an object in the bottom bar, then tap/click the tank to drop it.
- Drag on the liquid to stir it.
- Drag an object to grab and throw it.
- Tilt slider, left/right arrow keys, or the Motion button (phones with an accelerometer) tilt the tank.
- Gravity slider, Clear objects, Reset.

## How it works

- `engine.js`: a particle/grid (MPM-flavoured) liquid with five density layers. Rigid bodies sample the fluid at points across their shape to get buoyancy and drag, and push back on the fluid so waves and sloshing move the objects (and objects disturb the liquid).
- `index.html`: canvas rendering, controls and pointer handling.

The tank is built to fill the whole window: the grid size follows the screen, and the liquid is made of thousands of small particles (`TARGET_PARTICLES` in `index.html`). Each particle's mass is scaled by 2/ppc, so using more, lighter particles gives a finer liquid without changing its density or what floats.

No build step and no dependencies. Open `index.html` or serve the folder as a static site.

## Adding a new object

Add an entry to `SPECS` in `engine.js` (polygon outline plus either `dens` or `mass`, and a centre-of-mass offset), then add a button in `index.html` and a drawing branch in `drawBody`.
