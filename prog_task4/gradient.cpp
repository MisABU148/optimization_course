#include <iostream>
#include <iomanip> 

// f(x) = (-x)^2 + 4x + 1
// Function to maximize using the Gradient Ascent Method
float f(float x) {
    return -x * x + 4 * x + 1;
}

// f'(x) = -2x + 4
// The derivative is used to determine the gradient at any point x
float f_derivative(float x) {
    return -2 * x + 4;
}

// This function performs the Gradient Ascent algorithm to find the maximum of f(x).
//
// Parameters:
//
// x0: initial guess 
// alpha: learning rate 
// N: number of iterations
void gradient_ascent_method(float x0, float alpha, int N) {
    float x = x0;

    for (int i = 1; i <= N; i++) {
        float gradient = f_derivative(x);  
        x += alpha * gradient;             
    }

    std::cout << "Gradient Ascent Method:" << "\n";
    std::cout << "(xmax): " << std::fixed << std::setprecision(9) << x << ", f(xmax): " << f(x) << "\n";
}

int main() {
    // Input parameters
    float x0 = 0;           // initial guess 
    float alpha = 0.1;      // learning rate 
    int N = 100;            // number of iterations

    // Run Gradient Ascent with the given parameters
    gradient_ascent_method(x0, alpha, N);

    return 0; 
}
