#include <iostream>
#include <cmath>
#include <iomanip>

// Define the unimodal function
double f(double x) {
    return (x - 2) * (x - 2) + 3;
}

int main() {
    double a, b, epsilon;
    std::cout << "Enter interval [a, b]: ";
    std::cin >> a >> b;
    std::cout << "Enter tolerance epsilon: ";
    std::cin >> epsilon;
    
    // Golden ratio
    const double phi = (1 + std::sqrt(5.0)) / 2.0;

    // Start the Golden Section Search
    while ((b - a) > epsilon) {
        double x1 = b - (b - a) / phi;
        double x2 = a + (b - a) / phi;

        double f1 = f(x1);
        double f2 = f(x2);

        if (f1 >= f2) {
            // Minimum is in [x1, b]
            a = x1;
        } else {
            // Minimum is in [a, x2]
            b = x2;
        }
    }

    double xmin = (a + b) / 2.0;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Approximate xmin: " << xmin << std::endl;
    std::cout << "f(xmin): " << f(xmin) << std::endl;

    return 0;
}
