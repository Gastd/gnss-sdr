 /*!
 * \file beidou_sdr_signal_processing.cc
 * \brief This class implements various functions for BeiDou B1I signals
 * \author Giorgio Savastano, 2015. giorgio.savastano(at)uniroma1.it
 *
 * Detailed description of the file here if needed.
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 2010-2015  (see AUTHORS file for a list of contributors)
 *
 * GNSS-SDR is a software defined Global Navigation
 *          Satellite Systems receiver
 *
 * This file is part of GNSS-SDR.
 *
 * GNSS-SDR is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GNSS-SDR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNSS-SDR. If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include "beidou_b1i_signal_processing.h"
#include "BEIDOU_B1I.h"
#include <stdlib.h>
#include <cmath>

auto auxCeil = [](float x){ return static_cast<int>(static_cast<long>((x)+1)); };

auto mod = [](double a, double N){ return static_cast<int>(a - N*floor(a/N)); };

void beidou_b1i_code_gen_complex(std::complex<float>* _dest, signed int _prn, unsigned int _chip_shift)
{
    const unsigned int _code_length = BEIDOU_B1I_CODE_LENGTH_CHIPS;
    bool G1[_code_length];
    bool G2[_code_length];
    bool G1_register[11] = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
    bool G2_register[11] = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0};
    bool feedback1, feedback2;
    bool aux;
    unsigned int delay;
    unsigned int lcv, lcv2;

    short int prn_map[37][2] = {    1, 3,    // SV 1
                                    1, 4,    // SV 2
                                    1, 5,    // SV 3
                                    1, 6,    // SV 4
                                    1, 8,    // SV 5
                                    1, 9,    // SV 6
                                    1, 10,   // SV 7
                                    1, 11,   // SV 8
                                    2, 7,    // SV 9
                                    3, 4,    // SV 10
                                    3, 5,    // SV 11
                                    3, 6,    // SV 12
                                    3, 8,    // SV 13
                                    3, 9,    // SV 14
                                    3, 10,   // SV 15
                                    3, 11,   // SV 16
                                    4, 5,    // SV 17
                                    4, 6,    // SV 18
                                    4, 8,    // SV 19
                                    4, 9,    // SV 20
                                    4, 10,   // SV 21
                                    4, 11,   // SV 22
                                    5, 6,    // SV 23
                                    5, 8,    // SV 24
                                    5, 9,    // SV 25
                                    5, 10,   // SV 26
                                    5, 11,   // SV 27
                                    6, 8,    // SV 28
                                    6, 9,    // SV 29
                                    6, 10,   // SV 30
                                    6, 11,   // SV 31
                                    8, 9,    // SV 32
                                    8, 10,   // SV 33
                                    8, 11,   // SV 34
                                    9, 10,   // SV 35
                                    9, 11,   // SV 36
                                    10, 11}; // SV 37

    /* Generate G1 */
    for(lcv = 0; lcv < _code_length; lcv++)
        {
            // equal to the last value of the shift register
            G1[lcv] = G1_register[0];
            G2[lcv] = G2_register[prn_map[_prn -1][0] - 1] ^ G2_register[prn_map[_prn - 1][1] - 1];
            // computation of the G1 feedback
            feedback1 = (G1_register[0] + G1_register[1] + G1_register[2] + G1_register[3] + G1_register[4] + G1_register[10]) & 0x1;
            // computation of the G2 feedback
            feedback2 = (G2_register[0] + G2_register[2] + G2_register[3] + G2_register[6] + G2_register[7] + G2_register[8] + G2_register[9] + G2_register[10]) & 0x1;

            // shift
            for(lcv2 = 0; lcv2 < 10; lcv2++)
            {
                G1_register[lcv2] = G1_register[lcv2 + 1];
                G2_register[lcv2] = G2_register[lcv2 + 1];
            }
            // put feedback in position 1
            G1_register[10] = feedback1;
            G2_register[10] = feedback2;
        }

    /* Set the delay */
    // delay = _chip_shift;
    // delay %= _code_length;

    /* Generate PRN from G1 and G2 Registers */
    for (lcv = 0; lcv < _code_length; lcv++)
        {
            delay = (lcv + _chip_shift) % _code_length;
            aux = G1[delay] ^ G2[delay];
            if (aux == true)
                {
                    _dest[lcv] = std::complex<float>(1, 0);
                }
            else
                {
                    _dest[lcv] = std::complex<float>(-1, 0);
                }
            // delay++;
            // delay %= _code_length;
        }
}

/*
 *  Generates complex BeiDou B1I code for the desired SV ID and sampled to specific sampling frequency
 */
void beidou_b1i_code_gen_complex_sampled(std::complex<float>* _dest, unsigned int _prn, signed int _fs, unsigned int _chip_shift)
{
    // This function is based on the GNU software GPS for MATLAB in the Kay Borre book
    std::complex<float> _code[2046];
    signed int _samplesPerCode, _codeValueIndex;
    float _ts;
    float _tc;
    float aux;
    const signed int _codeFreqBasis = 2046000;  //Hz
    const signed int _codeLength = 2046;
    signed int _offset_prn, _offset_nh;
    double _phi_prn, _phi_nh;

    // const double _fs_in                 = static_cast<double>(_fs);
    // const signed int _codeLength        = static_cast<int>(BEIDOU_B1I_CODE_LENGTH_CHIPS);
    const signed int _codeDelayChips    = (_codeLength - _chip_shift) % _codeLength;
    const signed int _codeDelaySamples  = static_cast<int>(_codeDelayChips * (static_cast<double>(_fs) / BEIDOU_B1I_CODE_RATE_HZ));

    //--- Find number of samples per spreading code ----------------------------
    _samplesPerCode = static_cast<signed int>(static_cast<double>(_fs) / static_cast<double>(_codeFreqBasis / _codeLength));

    //generate B1I code 1 sample per chip
    beidou_b1i_code_gen_complex(_code, _prn, _chip_shift);

    for (signed int i = 0; i < _samplesPerCode; i++)
        {
            // Offset for the PRN codes in order to add the proper phase 
            _phi_prn = static_cast<double>(i - _codeDelaySamples) * (BEIDOU_B1I_CODE_RATE_HZ / static_cast<double>(_fs));
            _offset_prn = mod(_phi_prn, BEIDOU_B1I_CODE_LENGTH_CHIPS);

            // Offset for the NH code in order to add the proper phase
            _phi_nh = static_cast<double>(i - _codeDelaySamples) * (NH_BITS_RATE / static_cast<double>(_fs));
            _offset_nh = mod(_phi_nh, NH_BIT_DURATION);

            _dest[i] = _code[_offset_prn] * static_cast<float>(NH_CODE[_offset_nh]);
        }
}




