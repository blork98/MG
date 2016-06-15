#include "test.h"

#include <AutoLink.h>
#include <SparseMatrix.h>
#include <LinearIterativeSolver.h>
#include <GridOperators.h>
#include <MGTwoGridCorrection.h>

#include <vector>
#include <iostream>

void test_sparse()
{	
	std::vector<double> val = { 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0 };
	std::vector<unsigned int> colInd = { 0,3,0,1,3,0,2,3,4,2,3,4 };
	std::vector<unsigned int> rowPtr = { 0,2,5,9,11,12};
	size_t numCols = 5;
	size_t numRows = 5;

	SparseMatrix A(numRows, numCols, val, colInd, rowPtr);

	std::vector<double> sol(numRows, 0.0);
	std::vector<double> rhs = { 1.0,2.0,3.0,4.0,6.0 };


	//Test Diagonals
	A.get_diagonals(sol);
	std::cout << "Test Diagonals" << std::endl;
	for (auto it = sol.begin(); it != sol.end(); ++it)
		std::cout << *it << std::endl;

	//Test Vec Row Multiply
	std::cout << "Test Row Vec Multiply" << std::endl;
	for (size_t i = 0; i < numRows; ++i)
	{
		double res = A.vec_row_mult(i, rhs);
		std::cout << res << std::endl;
	}

	//Test Matrix Vec Mult
	std::cout << "Test Mat Vec mult" << std::endl;
	A.vec_multiply(rhs, sol);
	for (auto it = sol.begin(); it != sol.end(); ++it)
		std::cout << *it << std::endl;
};

void test_solver()
{
	std::vector<double> val = { 9.0,2.0,3.0,8.0,5.0,6.0,13.0,8.0,9.0,10.0,14.0,12.0 };
	std::vector<unsigned int> colInd = { 0,3,0,1,3,0,2,3,4,2,3,4 };
	std::vector<unsigned int> rowPtr = { 0,2,5,9,11,12 };
	size_t numCols = 5;
	size_t numRows = 5;

	SparseMatrix A(numRows, numCols, val, colInd, rowPtr);

	std::vector<double> sol(numRows, 1.0);
	std::vector<double> b = { 17.0,39.0,122.0,86.0,60.0 };

	std::shared_ptr<SparseMatrix> Mat = std::make_shared<SparseMatrix>(A);
	std::shared_ptr<std::vector<double>> rhs = std::make_shared<std::vector<double>>(b);

	double tol = 0.000001;
	unsigned int maxIter = 100;

	//Test Jacobi Solver
	JacobiSolver jacSolver(maxIter, tol);
	jacSolver.set_A(Mat);
	jacSolver.set_rhs(rhs);
	jacSolver.solve(sol);

	std::cout << "Test Jacobi Solver" << std::endl;
	for (auto it = sol.begin(); it != sol.end(); ++it)
		std::cout << *it << std::endl;

	//Test Gauss-Seidel Solver
	sol = std::vector<double>(numRows, 1.0);

	GaussSeidelSolver gsSolver(maxIter, tol);
	gsSolver.set_A(Mat);
	gsSolver.set_rhs(rhs);
	gsSolver.solve(sol);

	std::cout << "Test Gauss-Seidel Solver" << std::endl;
	for (auto it = sol.begin(); it != sol.end(); ++it)
		std::cout << *it << std::endl;

	//Test SOR
	sol = std::vector<double>(numRows, 1.0);
	double w = 1.2;
	bool isSym = false;
	SORSolver sorSolver(maxIter, tol, w, isSym);
	sorSolver.set_A(Mat);
	sorSolver.set_rhs(rhs);
	sorSolver.solve(sol);

	std::cout << "Test SOR Solver" << std::endl;
	for (auto it = sol.begin(); it != sol.end(); ++it)
		std::cout << *it << std::endl;

	//Test Symmetric SOR
	sol = std::vector<double>(numRows, 1.0);
	std::vector<double> val1 = { 2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0 };
	std::vector<unsigned int> colInd1 = { 0,1, 0,1,2, 1,2,3, 2,3,4, 3,4   };
	std::vector<unsigned int> rowPtr1 = { 0,2,5,8,11,13 };

	SparseMatrix Asym(numRows, numCols, val1, colInd1, rowPtr1);
	std::vector<double> b2 = { 0.0,0.0,0.0,0.0,6.0 };
	std::shared_ptr<SparseMatrix> Mat1 = std::make_shared<SparseMatrix>(Asym);
	std::shared_ptr<std::vector<double>> rhs1 = std::make_shared<std::vector<double>>(b2);

	isSym = true;
	w = 1.2;
	sorSolver.set_A(Mat1);
	sorSolver.set_rhs(rhs1);
	sorSolver.set_symmetric(isSym);
	sorSolver.solve(sol);
	sorSolver.set_w(w);

	std::cout << "Test Symmetric SOR Solver" << std::endl;
	for (auto it = sol.begin(); it != sol.end(); ++it)
		std::cout << *it << std::endl;
}

void test_grid_operators()
{
	size_t h = 8;
	size_t h2 = h/2;

	size_t fineGridSize = h + 1;
	size_t coarseGridSize = h2 + 1;
	std::vector<double> fineGrid = {0,1.7,2.5,3.5,4.2,3.7,1.2,4.5,1.8};
	std::vector<double> coarseGrid = {1.5,2.4,3.2,2.65,4.7 };
	std::vector<double> interpolatedFineGrid(fineGridSize, 0.0), restrictedCoarseGrid(coarseGridSize, 0.0);

	//Test Injection
	Injecction1D inject;
	inject.apply_operator(fineGrid, restrictedCoarseGrid);
	std::cout << "Test Restriction: Injection 1D" << std::endl;
	for (auto it = restrictedCoarseGrid.begin(); it != restrictedCoarseGrid.end(); ++it)
		std::cout << *it << std::endl;

	//Test Full Weighing
	restrictedCoarseGrid = std::vector<double>(coarseGridSize, 0.0);
	FullWeighing1D weighin;
	weighin.apply_operator(fineGrid, restrictedCoarseGrid);
	std::cout << "Test Restriction: Full Weighting 1D" << std::endl;
	for (auto it = restrictedCoarseGrid.begin(); it != restrictedCoarseGrid.end(); ++it)
		std::cout << *it << std::endl;

	//Test Linear Interpolation
	LinearInterpolation1D interp;
	interp.apply_operator(interpolatedFineGrid, coarseGrid);
	std::cout << "Test Interpolation: Linear Interpolation 1D" << std::endl;
	for (auto it = interpolatedFineGrid.begin(); it != interpolatedFineGrid.end(); ++it)
		std::cout << *it << std::endl;
}

void test_two_grid()
{
	//global vars
	int fineSweeps = 10000;
	double tolerance = 0.0001;
	double h = 1.0, scale = 1 / ((2 * h)*(2 * h));

	//Create coarse grid matrix
	std::vector<double> val1 = { 2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0 };
	std::vector<unsigned int> colInd1 = { 0,1, 0,1,2, 1,2,3, 2,3,4, 3,4 };
	std::vector<unsigned int> rowPtr1 = { 0,2,5,8,11,13 };

	size_t numRowsC = 5;
	SparseMatrix coarse(numRowsC, numRowsC, val1, colInd1, rowPtr1);
	coarse *= scale;
	std::shared_ptr<SparseMatrix> coarseA = std::make_shared<SparseMatrix>(coarse);

	//Create fine grid matrix
	std::vector<double> val2 = { 2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0,-1.0, 
		                        -1.0,2.0,-1.0, -1.0,2.0,-1.0, -1.0,2.0};
	std::vector<unsigned int> colInd2 = { 0,1, 0,1,2, 1,2,3, 2,3,4, 3,4,5, 4,5,6, 5,6,7, 6,7,8, 7,8 };
	std::vector<unsigned int> rowPtr2 = { 0,2,5,8,11,14,17,20,23,25 };

	size_t numRowsFine = 9;
	SparseMatrix fine(numRowsFine, numRowsFine, val2, colInd2, rowPtr2);
	std::shared_ptr<SparseMatrix> fineA = std::make_shared<SparseMatrix>(fine); 

	//Create Solvers
	std::shared_ptr<SparseIterativeLinearSolver> fineSolver(new JacobiSolver(fineSweeps, tolerance));
	std::shared_ptr<SparseLinearSolver> coarseSolver(new JacobiSolver(fineSweeps, tolerance));

	//Create rhs
	std::vector<double> fineRHS = {0.5,-1.0,2.0,0.0,-1.5,1.9,-1.1,1.5,-0.3};
	std::vector<double> coarseRHS = {-0.1,-0.6,3.3,-0.6,0.8};
	std::shared_ptr<std::vector<double>> fRHS(&fineRHS), cRHS(&coarseRHS);

	//solutions
	std::vector<double> fineSol(fineRHS.size(),0.0);
	std::vector<double> coarseSol(coarseRHS.size(), 0.0);

	/*
	//test solvers on coarse matrix
	coarseSolver->set_A(coarseA);
	coarseSolver->set_rhs(cRHS);
	coarseSolver->solve(coarseSol);
	std::cout << "Test Coarse Grid Matrix" << std::endl;
	for (auto it = coarseSol.begin(); it != coarseSol.end(); ++it)
		std::cout << *it << std::endl;
	std::cout << "Num Iters:" << std::dynamic_pointer_cast<SparseIterativeLinearSolver>(coarseSolver)->num_iters() << std::endl;

	//test solvers on fine matrix
	fineSolver->set_A(fineA);
	fineSolver->set_rhs(fRHS);
	fineSolver->solve(fineSol);
	std::cout << "Test Fine Grid Matrix" << std::endl;
	for (auto it = fineSol.begin(); it != fineSol.end(); ++it)
		std::cout << *it << std::endl;
	std::cout << "Num Iters:" << fineSolver->num_iters() << std::endl;
	*/

	//Declaret Interpolation, Restriction Operators
	std::shared_ptr<RestrictionOperator> R = std::make_shared<FullWeighing1D>(FullWeighing1D());
	std::shared_ptr<InterpolationOperator> I = std::make_shared<LinearInterpolation1D>(LinearInterpolation1D());

	//Create 2 Grid Multigrid
	TwoGridCorrection mg(fineSweeps, coarseSolver, fineSolver);
	mg.set_A_coarse(coarseA);
	mg.set_A_fine(fineA);
	mg.set_rhs(fRHS);
	mg.set_interpolation(I);
	mg.set_restriction(R);

	mg.solve(false, fineSol);
	std::cout << "Test Two Grid MG" << std::endl;
	for (auto it = fineSol.begin(); it != fineSol.end(); ++it)
		std::cout << *it << std::endl;

}
