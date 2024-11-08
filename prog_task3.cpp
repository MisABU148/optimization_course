#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

// Function to display the initial data in table format
void displayInitialData(const vector<int>& S, const vector<vector<int>>& C, const vector<int>& D) {
    int s = S.size(), d = D.size();

    cout << "Initial Data Table:" << endl;
    cout << setw(9) << "Supply" << " | ";
    for (int j = 0; j < d; j++) {
        cout << setw(5) << "D" << j + 1 << " ";
    }
    cout << endl;

    for (int i = 0; i < s; i++) {
        cout<< "S" << i + 1 << " = " << setw(4) << S[i] << " | ";
        for (int j = 0; j < d; j++) {
            cout << setw(6) << C[i][j] << " ";
        }
        cout << endl;
    }

    cout << setw(9) <<"Demands" << " | ";
    for (int j = 0; j < d; j++) {
        cout << setw(6) << D[j] << " ";
    }
    cout << endl << endl;
}

int Vogel_iteration(vector<int>& S, vector<vector<int>>& C, vector<int>& D, vector<vector<int>>& X0) {
    int s = S.size(), d = D.size();
    int max_dif = 0, min_row = -1, min_col = -1;

    // Calculate penalties for rows
    for (int i = 0; i < s; i++) {
        if (S[i] == 0) continue;  // Skip rows with zero supply
        int mn_1 = 1000, mn_2 = 1000;
        for (int j = 0; j < d; j++) {
            if (D[j] == 0) continue;  // Skip columns with zero demand
            if (C[i][j] < mn_1) {
                mn_2 = mn_1;
                mn_1 = C[i][j];
            } else if (C[i][j] < mn_2) {
                mn_2 = C[i][j];
            }
        }
        if (max_dif < abs(mn_2 - mn_1)) {
            max_dif = abs(mn_2 - mn_1);
            min_row = i;
        }
    }

    // Calculate penalties for columns
    for (int j = 0; j < d; j++) {
        if (D[j] == 0) continue;  // Skip columns with zero demand
        int mn_1 = 1000, mn_2 = 1000;
        for (int i = 0; i < s; i++) {
            if (S[i] == 0) continue;  // Skip rows with zero supply
            if (C[i][j] < mn_1) {
                mn_2 = mn_1;
                mn_1 = C[i][j];
            } else if (C[i][j] < mn_2) {
                mn_2 = C[i][j];
            }
        }
        if (max_dif < abs(mn_2 - mn_1)) {
            max_dif = abs(mn_2 - mn_1);
            min_col = j;
        }
    }

    // Find minimum cost in the selected row or column
    int min_weight = 1000;
    if (min_col == -1) {
        for (int j = 0; j < d; j++) {
            if (D[j] > 0 && C[min_row][j] < min_weight) {
                min_col = j;
                min_weight = C[min_row][j];
            }
        }
    } else {
        for (int i = 0; i < s; i++) {
            if (S[i] > 0 && C[i][min_col] < min_weight) {
                min_row = i;
                min_weight = C[i][min_col];
            }
        }
    }

    // Calculate allocation and update supply and demand
    int allocation = min(S[min_row], D[min_col]);
    X0[min_row][min_col] = allocation;  // Save the amount allocated

    int cost = allocation * C[min_row][min_col];
    S[min_row] -= allocation;
    D[min_col] -= allocation;

    return cost;
}

void Vogel(vector<int> S, vector<vector<int>> C, vector<int> D) {
    int total_cost = 0;
    vector<int> x;
    vector<vector<int>> X0(S.size(), vector<int>(D.size()));

    // Check for non-zero supply and demand left
    while (*max_element(S.begin(), S.end()) > 0 && *max_element(D.begin(), D.end()) > 0) {
        total_cost += Vogel_iteration(S, C, D, X0);
    }

    // Output results
    cout << endl << "Vogel's approximation Method:" << endl;
    for (int i = 0; i < S.size(); i++) {
        cout<< "X" << i + 1 << " | ";
        for (int j = 0; j < D.size(); j++) {
            cout << setw(6) << X0[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << "Total distribution cost: " << total_cost << endl;
}

int Russell_iteration(vector<int>& S, vector<vector<int>>& C, vector<int>& D, vector<vector<int>>& X0) {
    int s = S.size(), d = D.size();
    vector<int> u(s), v(d);

    // Compute u_i = max_j c_ij for each row i
    for (int i = 0; i < s; i++) {
        int max_c = 0;
        for (int j = 0; j < d; j++) {
            if (C[i][j] > max_c) {
                max_c = C[i][j];
            }
        }
        u[i] = max_c;
    }

    // Compute v_j = max_i c_ij for each column j
    for (int j = 0; j < d; j++) {
        int max_c = 0;
        for (int i = 0; i < s; i++) {
            if (C[i][j] > max_c) {
                max_c = C[i][j];
            }
        }
        v[j] = max_c;
    }

    // Compute modified costs c_ij' = c_ij - u_i - v_j
    int max_opportunity_cost = INT_MIN;
    int selected_i = -1, selected_j = -1;

    for (int i = 0; i < s; i++) {
        for (int j = 0; j < d; j++) {
            if (S[i] > 0 && D[j] > 0) {  // Only consider cells with non-zero supply and demand
                int modified_cost = C[i][j] - u[i] - v[j];
                if (modified_cost < 0) {
                    int opportunity_cost = -modified_cost;
                    if (opportunity_cost > max_opportunity_cost) {
                        max_opportunity_cost = opportunity_cost;
                        selected_i = i;
                        selected_j = j;
                    }
                }
            }
        }
    }

    // If no negative modified cost, select the cell with minimum c_ij
    if (selected_i == -1) {
        int min_c = INT_MAX;
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < d; j++) {
                if (S[i] > 0 && D[j] > 0 && C[i][j] < min_c) { // Only consider non-zero supply and demand
                    min_c = C[i][j];
                    selected_i = i;
                    selected_j = j;
                }
            }
        }
    }

    // Allocate as much as possible to cell (selected_i, selected_j)
    int amount = min(S[selected_i], D[selected_j]);
    int cost = amount * C[selected_i][selected_j];

    // Record allocation in the allocation matrix
    X0[selected_i][selected_j] = amount;

    // Adjust supplies and demands
    S[selected_i] -= amount;
    D[selected_j] -= amount;

    return cost;
}

void Russell(vector<int> S, vector<vector<int>> C, vector<int> D) {
    int total_cost = 0;
    vector<vector<int>> X0(S.size(), vector<int>(D.size(), 0)); // Initialize allocation matrix with zeros

    // Iterate until supplies and demands are exhausted
    while (*max_element(S.begin(), S.end()) > 0 && *max_element(D.begin(), D.end()) > 0) {
        total_cost += Russell_iteration(S, C, D, X0);
    }

    // Output the structured allocation matrix
    cout << "\nRussell's approximation Method:\n";

    for (int i = 0; i < X0.size(); i++) {
        cout << "X" << i + 1 << " | ";
        for (int j = 0; j < X0[i].size(); j++) {
            cout << setw(6) << X0[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\nTotal distribution cost: " << total_cost << "\n";
}

void northWest(vector<vector<int>> C, vector<int> S, vector<int> D) {
    int currentX = 0;
    int currentY = 0;
    double totalCost = 0;
    vector<vector<int>> X0(S.size(), vector<int>(D.size(), 0)); // Initialize allocation matrix with zeros

    cout << "Northwest Corner Method:\n";

    while (currentY < S.size() && currentX < D.size()) {
        // Determine the minimum value between the current row's supply and the current column's demand
        int minValue = min(S[currentY], D[currentX]);

        // Store the allocation in the allocation matrix
        X0[currentY][currentX] = minValue;

        // Increase the total cost by the product of the cost and the allocated amount
        totalCost += C[currentY][currentX] * minValue;

        // Subtract the allocated amount from the current supply and demand
        S[currentY] -= minValue;
        D[currentX] -= minValue;

        // Move to the next column or row if one of them is depleted
        if (S[currentY] == 0) {
            currentY++;  // Move to the next row
        }
        if (D[currentX] == 0) {
            currentX++;  // Move to the next column
        }
    }

    // Output the structured allocation matrix
    cout << "\nNorthwest Corner Allocation:\n";

    for (int i = 0; i < X0.size(); i++) {
        cout << "Ч" << i + 1 << " | ";
        for (int j = 0; j < X0[i].size(); j++) {
            cout << setw(6) << X0[i][j] << " ";
        }
        cout << endl;
    }

    cout << "\nTotal distribution cost: " << totalCost << "\n";
}

int main() {
    int s = 3, d = 4;
    int c, sum_s = 0, sum_d = 0;
    vector<int> S;
    for (int i = 0; i < s; i++) {
        cin >> c;
        sum_s += c;
        S.push_back(c);
    }

    vector<vector<int>> C(s, vector<int>(d));
    for (int i = 0; i < s; i++) {
        for (int j = 0; j < d; j++) {
            cin >> C[i][j];
        }
    }

    vector<int> D;
    for (int i = 0; i < d; i++) {
        cin >> c;
        sum_d += c;
        D.push_back(c);
    }

    if (sum_s != sum_d) {
        cout << "The problem is not balanced!";
        return 0;
    }

    // Display the initial data table
    displayInitialData(S, C, D);

    northWest(C, S, D);

    //Volgel`s approximation
    Vogel(S, C, D);

    //Russell`s approximation
    Russell(S, C, D);

    return 0;
}
