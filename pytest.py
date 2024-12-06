#!/usr/bin/env python
from ctypes import *

from scipy.signal import windows


class Test(object):
    EPSILON = 1e-14
    WIN_SIZE = 512
    dll = CDLL('./libwf.dll')
    test_arr = (c_double * WIN_SIZE)()

    def __init__(self, dll_wf_name: str, scipy_wf_name: str, has_alpha=False, alpha=0.0):
        self.dll_wf_name = dll_wf_name
        self.scipy_wf_name = scipy_wf_name
        self.has_alpha = has_alpha
        self.alpha = alpha

    def test(self) -> None:
        res = True

        print(f'test for "{self.dll_wf_name}({
              self.alpha if self.has_alpha else ''})"')
        # функции для получения массивов
        wf = getattr(self.dll, self.dll_wf_name)
        scipy_wf = getattr(windows, self.scipy_wf_name)
        # создание массивов
        if not self.has_alpha:
            wf(self.test_arr, c_size_t(self.WIN_SIZE))
            test_win = scipy_wf(self.WIN_SIZE)
        else:
            wf(self.test_arr, c_size_t(self.WIN_SIZE), c_double(self.alpha))
            test_win = scipy_wf(self.WIN_SIZE, self.alpha)
        # сравнение массивов
        for i in range(self.WIN_SIZE):
            xn = float(self.test_arr[i])
            ref = float(test_win[i])

            if abs(xn - ref) > self.EPSILON:
                res = False
                print(f'[{i:03}]: {xn} != {ref}')

        return res


if __name__ == '__main__':
    TESTS = [
        Test('wf_rect', 'boxcar'),
        Test('wf_triang', 'triang'),
        Test('wf_bartlett', 'bartlett'),

        Test('wf_hann', 'hann'),
        Test('wf_hamming', 'hamming'),
        Test('wf_general_hamming', 'general_hamming', True, 0.14),
        Test('wf_general_hamming', 'general_hamming', True, 0.88),
        Test('wf_blackman', 'blackman'),
        Test('wf_nuttall', 'nuttall'),
        Test('wf_blackmanharris', 'blackmanharris'),
        Test('wf_flattop', 'flattop'),

        Test('wf_tukey', 'tukey', True, 0.3),
        Test('wf_tukey', 'tukey', True, 0.5),
        Test('wf_tukey', 'tukey', True, 0.8),
    ]

    for t in TESTS:
        t.test()
