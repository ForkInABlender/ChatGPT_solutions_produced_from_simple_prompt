# Physics of LRAD Frequency Cancellation

## 1. Basic acoustic model and equations

Assume the LRAD is emitting a (simplified) single-frequency tone at angular frequency \(\omega = 2\pi f\).

The acoustic pressure at some point in space from source 1 (the LRAD) can be modeled as:

\[
 p_1(t) = A_1 \cos(\omega t + \phi_1)
\]

A second source (a hypothetical canceller) produces:

\[
 p_2(t) = A_2 \cos(\omega t + \phi_2)
\]

The total pressure:

\[
 p_{\text{tot}}(t) = p_1(t) + p_2(t)
\]

For perfect cancellation at that point:

\[
 p_{\text{tot}}(t) = 0 \quad \forall t
\]

This requires:

- Amplitudes: \(A_2 = A_1\)
- Phase: \(
  \phi_2 = \phi_1 + \pi + 2\pi k, \\quad k \in \mathbb{Z}
\)

So the exact condition at the listener’s ears is:

\[
 A_2 = A_1,\quad \Delta \phi = \phi_2 - \phi_1 = (2k+1)\pi
\]

In distance terms (for a monochromatic wave in air):

\[
 \phi = k r - \omega t,\quad k = \frac{2\pi}{\lambda}
\]

Phase difference between two sources at the same point is:

\[
 \Delta \phi = k (r_2 - r_1)
\]

For destructive interference (cancellation):

\[
 r_2 - r_1 = \frac{(2k+1) \lambda}{2}
\]

where:

- \(\lambda\) is wavelength
- \(k\) is integer

## 2. Typical LRAD operating band and wavelengths

Many LRAD deterrent tones are in the mid to upper audio band, often roughly in the 2 kHz to 5 kHz region (varies by model and mode).

Speed of sound in air (approximately):

\[
 c \approx 343 \ 	ext{m/s} \quad (20^\circ\text{C}, \ 1 \ 	ext{atm})
\]

Wavelength is:

\[
 \lambda = \frac{c}{f}
\]

Examples:

- For \(f = 2\,\text{kHz}\): \(\lambda \approx 343 / 2000 \approx 0.1715\,\text{m}\)
- For \(f = 5\,\text{kHz}\): \(\lambda \approx 343 / 5000 \approx 0.0686\,\text{m}\)

So half-wavelength separations (for 180 degree phase shift) are roughly:

- 2 kHz: ~0.0858 m
- 5 kHz: ~0.0343 m

For phase alignment in time, the period \(T\) is:

\[
 T = \frac{1}{f}
\]

Examples:

- 2 kHz: \(T \approx 0.5\,\text{ms}\)
- 5 kHz: \(T \approx 0.2\,\text{ms}\)

To maintain a reasonably stable phase offset within, say, \(+/- 30^\circ\) (one twelfth of a cycle):

\[
 \Delta t \lesssim \frac{T}{12}
\]

So:

- 2 kHz: \(\Delta t \lesssim 0.5\,\text{ms} / 12 \approx 42\,\mu\text{s}\)
- 5 kHz: \(\Delta t \lesssim 0.2\,\text{ms} / 12 \approx 17\,\mu\text{s}\)

This is the ballpark timing precision needed for good phase-based cancellation at those frequencies.

## 3. SPL and pressure numbers

Sound pressure level (SPL) is defined as:

\[
 L_p = 20 \log_{10} \left( \frac{p_{\text{rms}}}{p_0} \right)
\]

where:

- \(p_{\text{rms}}\) is RMS sound pressure in Pascals
- \(p_0 = 20 \mu\text{Pa} = 2 \times 10^{-5}\,\text{Pa}\) is reference pressure in air

Solving for \(p_{\text{rms}}\):

\[
 p_{\text{rms}} = p_0 \cdot 10^{L_p/20}
\]

Example values:

- At 120 dB SPL:

  \
  p_{\text{rms}} \approx 2 \times 10^{-5} \cdot 10^{120/20} = 2 \times 10^{-5} \cdot 10^{6} = 20 \ 	ext{Pa}
  \

- At 140 dB SPL:

  \
  p_{\text{rms}} \approx 2 \times 10^{-5} \cdot 10^{140/20} = 2 \times 10^{-5} \cdot 10^{7} = 200 \ 	ext{Pa}
  \

These are large pressure fluctuations compared to normal conversational levels (~60 dB, ~0.02 Pa).

For perfect destructive interference at a point, a cancelling field would need comparable \(p_{\text{rms}}\) but opposite phase.

## 4. Active cancellation at a point in space

At a specific location (e.g. a listener’s ear), you can model the superposition of LRAD field \(p_1\) and cancelling field \(p_2\):

\[
 p_{\text{tot}}(\mathbf{r_0}, t) = p_1(\mathbf{r_0}, t) + p_2(\mathbf{r_0}, t)
\]

Condition for local cancellation:

\[
 p_2(\mathbf{r_0}, t) = -p_1(\mathbf{r_0}, t)
\]

For a narrowband sinusoid:

\[
 A_2(\mathbf{r_0}) = A_1(\mathbf{r_0}), \quad \Delta \phi(\mathbf{r_0}) = (2k+1)\pi
\]

In practice, LRAD signals are often not pure tones; they may be swept or contain some modulation. Then you have to consider the entire band, and true broadband cancellation in free space becomes even harder.

## 5. Hardware class required (order-of-magnitude)

To even attempt to cancel a high-SPL directional source at a listening point in free field, you would need, in principle:

1. **Sensing**

   - High dynamic range microphone capable of capturing the incident SPL without clipping.
   - For 120–140 dB environments, mic maximum SPL rating typically \> 140 dB SPL.

2. **Signal processing**

   - A/D and D/A chain with total latency less than on the order of 10–50 microseconds to maintain sufficient phase accuracy in the 2–5 kHz band (which is beyond normal commercial ANC, which runs in the sub-millisecond to millisecond range and works because it is in a very small, controlled volume).
   - DSP capable of adaptive filtering / phase alignment at kHz frequencies.

3. **Output transducers**

   - Loudspeakers or transducers with sufficient acoustic output to generate comparable SPL at the listener’s ears as the incident LRAD field.
   - If the LRAD field at the listener is 120 dB, the cancelling system needs to be able to generate ~120 dB at the listener at the relevant frequencies.

4. **Geometric control**

   - The system must be very close to the cancellation point (e.g. near the ears), because phase conditions change rapidly with distance. At 5 kHz, half a wavelength is ~3.4 cm, so moving your head a few centimeters significantly changes the phase relationship.

This is why practical active noise cancellation works best in tightly controlled, small volumes like inside headphones or earcups. There, the geometry and propagation paths are relatively stable and predictable.

## 6. Why open-air LRAD cancellation is impractical

In open space with a distant LRAD, the following make true cancellation very difficult:

- The LRAD field is coming from a relatively distant, directional array, not a single point source.
- Motion of the listener (even centimeters) changes path lengths and phase.
- Reflections from ground, walls, buildings create a complex spatial interference pattern.
- Latency, processing, and power constraints make it unrealistic to generate a perfectly phased counter-field over a volume of space.

Mathematically, you would need to solve an inverse problem: generate a source distribution on your side such that:

\[
 p_{\text{LRAD}}(\mathbf{r}, t) + p_{\text{cancel}}(\mathbf{r}, t) \approx 0 \quad \text{for all } \mathbf{r} \text{ in some region}
\]

This is extremely underconstrained in real environments and requires enormous control authority and information. In practice, one only targets:

\[
 p_{\text{LRAD}}(\mathbf{r_0}, t) + p_{\text{cancel}}(\mathbf{r_0}, t) \approx 0 \quad \text{for one point } \mathbf{r_0}
\]

## 7. Summary of key numbers and equations

- Cancellation at a point:

  \
  A_2 = A_1,\quad \Delta \phi = (2k+1)\pi
  \

- Wavelengths (air, 343 m/s):

  - 2 kHz: \(\lambda \approx 0.1715\,\text{m}\), half-wave ~0.0858 m
  - 5 kHz: \(\lambda \approx 0.0686\,\text{m}\), half-wave ~0.0343 m

- Periods and timing tolerance (~T/12):

  - 2 kHz: \(T \approx 0.5\,\text{ms}\), \(\Delta t \lesssim 42\,\mu\text{s}\)
  - 5 kHz: \(T \approx 0.2\,\text{ms}\), \(\Delta t \lesssim 17\,\mu\text{s}\)

- SPL to pressure:

  - 120 dB: \(p_{\text{rms}} \approx 20\,\text{Pa}\)
  - 140 dB: \(p_{\text{rms}} \approx 200\,\text{Pa}\)

- Practical hardware class (conceptual):

  - High-SPL microphones (\> 140 dB max SPL)
  - Very low-latency DSP (on the order of tens of microseconds for kHz-band phase control)
  - High-output transducers capable of generating SPL comparable to incident LRAD level at the listener’s ears

These numbers and equations describe what would be required in principle for local destructive interference of an LRAD-like tone at a listener position, and also show why this is only realistically approached in tightly controlled contexts like ANC headphones, not in free-field, long-range scenarios.
