
//   HyKiCT - 1D Lagrangian Radiation-Hydrodyanmics code for testing Coupling with Kinetic codes.
//   Copyright (C) 2020- Abetharan Antony
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <https://www.gnu.org/licenses/>. 
#ifndef FASTMATHS_HPP
#define FASTMATHS_HPP
 
inline float InvSqrt(float number)
{
    union {
            float f;
            uint32_t i;
        } conv;

        float x2;
        const float threehalfs = 1.5F;

        x2 = number * 0.5F;
        conv.f  = number;
        conv.i  = 0x5f3759df - ( conv.i >> 1 );
        conv.f  = conv.f * ( threehalfs - ( x2 * conv.f * conv.f ) );
        return conv.f;
}

#endif

