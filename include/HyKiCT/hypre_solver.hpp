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
#ifndef HYPRESOLVER_HPP
#define HYPRESOLVER_HPP
#include <numeric>
#include <algorithm>
#include "GridData.h"
#include "FixedData.h"
#include "PlasmaParameters.h"
#include "HYPRE.h"
#include "HYPRE_struct_ls.h
// #include "Switches.hpp"

class HypreSolver
{
    private:
            HYPRE_StructGrid     grid;
            HYPRE_StructStencil  stencil;
            HYPRE_StructMatrix   A;
            HYPRE_StructVector   b;
            HYPRE_StructVector   x;
            HYPRE_StructSolver   solver;
            HYPRE_StructSolver   precond;
    public:
        HypreSolver()
        {
            HYPRE_Init();
            HYPRE_StructGridCreate(MPI_COMM_WORLD, 1, &grid);
            
        }

};
#endif