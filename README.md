# BadFFT
Wrapper file for [idutny's fft.js](https://github.com/indutny/fft.js/)

Adds a function ``` badfft( data ) ``` to the global scope, which takes an array input and returns an array where each item is the log of the magnitude of that term in the fourier transform. The output is shifted so that the lowest frequencies appear in the middle, to create a centered spectrum.
