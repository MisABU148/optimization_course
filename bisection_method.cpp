#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// Define the function f(x)
double f(double x) {
    return x*x*x - 6*x*x + 11*x - 6; // f(x) = x^3 - 6x^2 + 11x - 6
}

int main() {
    double a, b, epsilon;
    
    // User input
    cout << "Enter the initial interval [a, b]: " << endl;
    cin >> a >> b;
    cout << "Enter the tolerance (epsilon): " << endl;
    cin >> epsilon;
    
    // Check if f(a)*f(b) < 0 to ensure a root lies in [a,b]
    if (f(a) * f(b) >= 0) {
        cout << "The function does not have opposite signs at the interval endpoints." << endl;
        cout << "Please choose another interval." << endl;
        return 1; // Exit the program with an error code
    }

    double c;
    int iteration = 0;
    
    cout << fixed << setprecision(10);
    cout << "Iter\t   a\t\t   b\t\t   c\t\t   f(c)" << endl;
    
    while (true) {
        iteration++;
        c = (a + b) / 2.0;
        double fc = f(c);

        // Print current iteration details
        cout << iteration << "\t" << a << "\t" << b << "\t" << c << "\t" << fc << endl;

        // Check if |f(c)| < epsilon
        if (fabs(fc) < epsilon) {
            cout << "Root found at c = " << c << " with f(c) = " << fc << endl;
            break;
        }

        // Decide which half of the interval to use for the next iteration
        if (f(a) * fc < 0) {
            b = c; 
        } else {
            a = c; 
        }
    }

    return 0;
}
