#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

// Helper function to compute the transpose
vector<vector<double>> transpose(const vector<vector<double>> &matrix)
{
    int rows = matrix.size();
    int cols = matrix[0].size();

    // Create a matrix to hold the transpose
    vector<vector<double>> transposedMatrix(cols, vector<double>(rows, 0.0));

    // Compute the transpose
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            transposedMatrix[j][i] = matrix[i][j];
        }
    }
    return transposedMatrix;
}

// Helper function to calculate the inverse of a matrix
vector<vector<double>> inverse(const vector<vector<double>> &matrix) {
    int n = matrix.size();

    // Create an augmented matrix [A | I] where I is the identity matrix
    vector<vector<double>> augmentedMatrix(n, vector<double>(2 * n, 0.0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            augmentedMatrix[i][j] = matrix[i][j]; // Copy original matrix
        }
        augmentedMatrix[i][n + i] = 1.0; // Identity matrix
    }

    // Perform Gaussian elimination with partial pivoting
    for (int i = 0; i < n; ++i)
    {
        // Find the pivot row
        int pivotRow = i;
        for (int j = i + 1; j < n; ++j)
        {
            if (fabs(augmentedMatrix[j][i]) > fabs(augmentedMatrix[pivotRow][i]))
            {
                pivotRow = j;
            }
        }

        // Swap rows if needed
        if (pivotRow != i)
        {
            swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);
        }

        // Normalize the pivot row
        double pivotValue = augmentedMatrix[i][i];
        if (fabs(pivotValue) < 1e-10)
        {
            throw runtime_error("Matrix is singular and cannot be inverted.");
        }
        for (int j = 0; j < 2 * n; ++j)
        {
            augmentedMatrix[i][j] /= pivotValue;
        }

        // Eliminate other rows
        for (int j = 0; j < n; ++j)
        {
            if (i != j)
            {
                double factor = augmentedMatrix[j][i];
                for (int k = 0; k < 2 * n; ++k)
                {
                    augmentedMatrix[j][k] -= factor * augmentedMatrix[i][k];
                }
            }
        }
    }

    // Extract the inverse matrix from the augmented matrix
    vector<vector<double>> inverseMatrix(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            inverseMatrix[i][j] = augmentedMatrix[i][n + j];
        }
    }

    return inverseMatrix;
}

// Helper function to multiply matrix by scalar
vector<vector<double>> multiplyMatrixByScalar(vector<vector<double>> matrix,
                                              double scalar) {
    // Get the dimensions of the matrix
    int rows = matrix.size();
    int cols = matrix[0].size();

    // Create a result matrix with the same dimensions
    vector<vector<double>> result(rows, vector<double>(cols));

    // Multiply each element by the scalar
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result[i][j] = matrix[i][j] * scalar;
        }
    }

    // Return the result matrix
    return result;
}

// Helper function to multiply two matrices A and B
vector<vector<double>> multiplyMatrices(vector<vector<double>> A,
                                        vector<vector<double>> B) {
    int rowsA = A.size();	 // Number of rows in matrix A
    int colsA = A[0].size(); // Number of columns in matrix A
    int rowsB = B.size();	 // Number of rows in matrix B
    int colsB = B[0].size(); // Number of columns in matrix B

    // Check if the matrices can be multiplied
    if (colsA != rowsB)
    {
        throw invalid_argument(
                "Number of columns in A must equal number of rows in B for "
                "multiplication.");
    }

    // Resultant matrix will have rowsA rows and colsB columns
    vector<vector<double>> result(rowsA, vector<double>(colsB, 0.0));

    // Multiply A and B
    for (int i = 0; i < rowsA; ++i)
    {
        for (int j = 0; j < colsB; ++j)
        {
            for (int k = 0; k < colsA; ++k)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

// Interior Point Algorithm function
void interiorPointMethod(vector<vector<double>> A, vector<double> b,
                         vector<double> c, double alpha)
{
    int n = c.size(); // Number of decision variables
    int m = b.size(); // Number of constraints

    for (int i = 0; i < m; ++i)
    {
        c.push_back(0);
    }

    // Initial guess for x, starting with positive values to be in the interior
    vector<double> x(n, 1.0); // Initialize with a small value like 1.0

    // Introduce artificial variables s to make the system feasible
    vector<double> s(m, 0.0);

    // Ensure initial feasibility
    for (int i = 0; i < m; ++i)
    {
        double sum = 0;
        for (int j = 0; j < n; ++j)
        {
            sum += A[i][j] * x[j];
        }

        // If constraint is violated, we add slack
        if (sum < b[i])
        {
            s[i] = b[i] - sum;
        }
        else
        {
            s[i] = 0.0; // No infeasibility, no slack needed
        }
    }

    x.insert(x.end(), s.begin(), s.end());

    double optimum;

    int maxIterations = 20;

    vector<vector<double>> A_initial = A;

    for (int iter = 0; iter < maxIterations; ++iter) {
//        cout << endl
//             << "ITERATION #" << iter << endl;

        // Compose D matrix
        vector<vector<double>> D;
        for (int i = 0; i < n + m; ++i)
        {
            D.push_back(vector<double>(n + m, 0));
        }

        for (int i = 0; i < n + m; ++i)
        {
            D[i][i] = x[i];
        }

        // Calculate A~ = AD
        A = multiplyMatrices(A_initial, D);

        // Calculate P = I - Transpose(A~) * Inverse(A~ * Transpose(A~)) * A~
        vector<vector<double>> I;
        for (int i = 0; i < n + m; ++i)
        {
            I.push_back(vector<double>(n + m, 0));
        }
        for (int i = 0; i < n + m; ++i)
        {
            I[i][i] = 1;
        }
        vector<vector<double>> A_trans = transpose(A);

        // Transpose(A~) * Inverse(A~ * Transpose(A~)) * A~
        vector<vector<double>> subtrahend = multiplyMatrices(
                multiplyMatrices(A_trans,
                                 inverse(multiplyMatrices(A, A_trans))
                ),
                A);

        vector<vector<double>> P;
        for (int i = 0; i < n + m; ++i)
        {
            P.push_back(vector<double>(n + m, 0));
        }

        for (int i = 0; i < n + m; ++i)
        {
            for (int j = 0; j < n + m; ++j)
            {
                P[i][j] = I[i][j] - subtrahend[i][j];
            }
        }

        // Calculate Cp = P * (Dc)
        vector<vector<double>> cp =
                multiplyMatrices(P, multiplyMatrices(D, transpose({c})));

        // Calculate nu - absolute value of the largest negative cp[i]
        cp = transpose(cp);
        double min = *min_element(cp[0].begin(), cp[0].end());
        double v = abs(min);

        cp = transpose(cp);

        // Calculate X~
        vector<double> ones;
        for (int i = 0; i < n + m; ++i)
        {
            ones.push_back(1);
        }

        vector<vector<double>> cp_with_coeff =
                multiplyMatrixByScalar(cp, (alpha / v));

        vector<vector<double>> X_tilde;
        for (int i = 0; i < n + m; ++i)
        {
            X_tilde.push_back({});
            X_tilde[i].push_back(1 + cp_with_coeff[i][0]);
        }

        // Calculate new X = D * X~
        x = transpose(multiplyMatrices(D, X_tilde))[0];

        bool inequlities_are_true = true;

        for (int i = 0; i < m; ++i)
        {
            double r = 0;
            for (int j = 0; j < n; ++j)
            {
                r += A_initial[i][j] * x[j];
//                cout << A_initial[i][j] << " * " << x[j] << "  +  ";
            }
//            cout << endl;
            inequlities_are_true = r < b[i];
//            cout << r << " < " << b[i] << " -> " << is_less << endl;
        }
        if (!inequlities_are_true)
        {
            cout << "Constraints are not satisfied, the end of iterations" << endl;
            break;
        }

        // Calculate and print optimum
        optimum = 0;
//        cout << "new X: ";
        for (int i = 0; i < n; ++i)
        {
//            cout << x[i] << " ";
            optimum += x[i] * c[i];
        }
//        cout << endl
//             << "OPTIMUM = " << optimum << endl;
    }

    cout << endl << "Answer for Interior-Point algorithm when a = " << alpha << " :" << endl;
    for (int i = 0; i < n; ++i)
    {
        cout << x[i] << " ";
    }
    cout << endl << "Optimum = " << optimum << endl;
}

bool is_appl = 1;

bool Simplex_method(vector<vector<double>> A, vector<double> b,
                    vector<double> c) {
    int n = A.size(), m = A[0].size();
    double A1[4][8] {};
    double mn = 0;
    int enter = 0, leave=0;
    for (int i = 1; i<m-1; i++) {
        if (A[0][i]<mn) {
            mn = A[0][i];
            enter = i;
        }
    }
    if (mn == 0)
        return 1;

    mn = 1000000;
    for (int i = 1; i<n; i++) {
        if (A[i][enter]>0 && A[i][m-1]/A[i][enter] < mn) {
            mn = A[i][m-1] / A[i][enter];
            leave = A[i][0];
        }
    }
    if (mn == 1000000) {
        is_appl = 0;
        return 1;
    }

    for (int i = 1; i<n; i++){
        if (A[i][0] == leave) {
            A[i][0] = enter;
            leave = i;
            break;
        }
    }



    for (int i = 0; i<n; i++) {
        if (A[i][0] == enter) {
            A1[i][0] = enter;
            for (int j = 1; j<m;j++) {
                if (A[i][enter] != 0)
                    A1[i][j] = A[i][j] / A[i][enter];
                else
                    A1[i][j] = 0;
            }
            continue;
        }
        A1[i][0] = A[i][0];
        for (int j = 1; j<m; j++) {
            if (j == enter) {
                A1[i][j] = (i == enter)?1:0;
            } else {
                if (A[i][enter] != 0)
                    A1[i][j] = A[i][j] - (A[i][enter] * A[leave][j] / A[leave][enter]);
                else
                    A1[i][j] = 0;
            }
        }
    }

    for (int i = 0; i < n;i++) {
        for (int j = 0; j < m;j++) {
            A[i][j] = A1[i][j];
        }
    }

    return 0;
}

int main()
{
    int n = 3, m = 6;
    double num;
    vector<double> c;
    for (int i = 0; i < 3;i++) {
        cin >> num;
        c.push_back(num);
    }
    vector<vector<double>> A;
    for (int i = 0; i < n;i++) {
        A.push_back({});
        for (int j = 0; j < m;j++) {
            cin >> num;
            A[i].push_back(num);
        }
    }
    vector<double> b;
    for (int i = 0; i < 3;i++) {
        cin >> num;
        b.push_back(num);
    }
    double eps;
    cin >> eps;
    double alpha1 = 0.5;
    double alpha2 = 0.9;

    interiorPointMethod(A, b, c, alpha1);
    interiorPointMethod(A, b, c, alpha2);

    while (!Simplex_method(A, b, c)) {}

    if (is_appl) {
        cout << endl << "Answer for Simplex method :" << endl;
        for (int i = 1; i < n;i++) {
            int in_ans = 0;
            for (int j = 1; j < n;j++) {
                if (A[j][0] == i) {
                    cout << A[j][m-1] << ' ';
                    in_ans = 1;
                    break;
                }
            }
            if(!in_ans)
                cout << 0 << ' ';
        }
        cout << endl <<
            "Optimum = " << A[0][m-1];
    } else {
        cout << "The method is not applicable!";
    }

    return 0;
}