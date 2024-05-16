from __future__ import annotations
import scikit_build_example as m
import numpy as np
import wave
import struct

def read_wave(filename: str) -> list[int]:
    with wave.open(filename) as f:
        data = f.readframes(f.getnframes())
        samples = struct.unpack('{n}h'.format(n=f.getnframes()*f.getnchannels()), data)
    return samples

audio = read_wave('file.wav')
dft_audio = []

m.visualization(audio)

dft = m.dft(audio)
dft_modulo = m.modulo(dft)

m.visualization(dft_modulo)

inverse = m.inverse_transform(dft)

m.visualization(inverse)