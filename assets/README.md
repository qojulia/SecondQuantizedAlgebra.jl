# SecondQuantizedAlgebra.jl logo

`logo.svg` is the canonical package logo. `docs/src/assets/logo.svg` is an
identical copy. `icon.svg` is the square favicon form, hinted to a 16-pixel grid.

The mark is the commutator `[ô,ô]` drawn as an owl:

- the brackets are plain square brackets drawn as constant-width strokes with
  round caps and joins, tall enough to enclose the hats;
- the operators are exact filled circles (the eyes), each with a circumflex hat
  in its own color and a small white catchlight at the same offset;
- the comma, a smooth 16-segment cubic path, doubles as the beak;
- the colors are Julia purple `#9558B2`, green `#389826`, and red `#CB3C33`,
  with neutral gray `#888888` brackets.

All logo assets are standalone SVGs: no fonts or raster images are embedded.
The canonical logo is used by Documenter as the sidebar logo and landing-page
hero. `docs/src/assets/favicon.ico` (16, 32, 48 and 64 px) is rendered from
`icon.svg`, since the wide logo does not survive at 16 px.
