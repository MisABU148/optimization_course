#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

int Vogel_iteration(vector<int>& S, vector<vector<int>>& C, vector<int>& D, vector<int>& X) {
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
    X.push_back(min_row + 1);
    X.push_back(min_col + 1);
    X.push_back(allocation);  // Save the amount allocated

    int cost = allocation * C[min_row][min_col];
    S[min_row] -= allocation;
    D[min_col] -= allocation;

    return cost;
}

void Vogel(vector<int> S, vector<vector<int>> C, vector<int> D) {
    int total_cost = 0;
    vector<int> x;

    // Check for non-zero supply and demand left
    while (*max_element(S.begin(), S.end()) > 0 && *max_element(D.begin(), D.end()) > 0) {
        total_cost += Vogel_iteration(S, C, D, x);
    }

    // Output results
    cout << "Vogel's Approximation Method:" << endl;
    for (int i = 0; i < x.size(); i += 3) {
        cout << "x_" << x[i] << x[i+1] << " = " << x[i+2] << ' ';
    }
    cout << endl << "Total distribution cost: " << total_cost << endl;
}

// Function to display the initial data in table format
void displayInitialData(const vector<int>& S, const vector<vector<int>>& C, const vector<int>& D) {
    int s = S.size(), d = D.size();

    cout << "Initial Data Table:" << endl;
    cout << setw(9) << "Supply" << " | ";
    for (int j = 0; j < d; j++) {
        cout << setw(4) << "D" << j + 1 << " ";
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

int Russell_iteration(vector<int>& S, vector<vector<int>>& C, vector<int>& D) {
    int s = S.size(), d = D.size();
    vector<int> u(s), v(d);

    // Compute u_i = min_j c_ij for each row i
    for (int i = 0; i < s; i++) {
        int max_c = 0;
        for (int j = 0; j < d; j++) {
            if (C[i][j] > max_c) {
                max_c = C[i][j];
            }
        }
        u[i] = max_c;
//        cerr << "U[" << i << "] = " << u[i] << "\n";
    }

    // Compute v_j = min_i c_ij for each column j
    for (int j = 0; j < d; j++) {
        int max_c = 0;
        for (int i = 0; i < s; i++) {
            if (C[i][j] > max_c) {
                max_c = C[i][j];
            }
        }
        v[j] = max_c;
//        cerr << "V[" << j << "] = " << v[j] << "\n";
    }

    // Compute modified costs c_ij' = c_ij - u_i - v_j
    int max_opportunity_cost = INT_MIN;
    int selected_i = -1, selected_j = -1;

    for (int i = 0; i < s; i++) {
        for (int j = 0; j < d; j++) {
            int modified_cost = C[i][j] - u[i] - v[j];
//            cerr << modified_cost << " ";
            if (modified_cost < 0) {
                int opportunity_cost = -modified_cost;
                if (opportunity_cost > max_opportunity_cost) {
                    max_opportunity_cost = opportunity_cost;
                    selected_i = i;
                    selected_j = j;
                }
            }
        }
//        cerr << "\n";
    }
//    cerr << "--\n";

    // If no negative modified cost, select the cell with minimum c_ij
    if (selected_i == -1) {
        int min_c = INT_MAX;
        for (int i = 0; i < s; i++) {
            for (int j = 0; j < d; j++) {
                if (C[i][j] < min_c) {
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

    // Adjust supplies and demands
    S[selected_i] -= amount;
    D[selected_j] -= amount;

    // If supply exhausted, remove row
    if (S[selected_i] == 0) {
        S.erase(S.begin() + selected_i);
        C.erase(C.begin() + selected_i);
    }

    // If demand exhausted, remove column
    if (D[selected_j] == 0) {
        D.erase(D.begin() + selected_j);
        for (int i = 0; i < S.size(); i++) {
            C[i].erase(C[i].begin() + selected_j);
        }
    }

    return cost;
}

int Russell(vector<int> S, vector<vector<int>> C, vector<int> D) {
    int total_cost = 0;
    while (!S.empty() && !D.empty()) {
        total_cost += Russell_iteration(S, C, D);
    }
    return total_cost;
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

    //Volgel`s approximation
    Vogel(S, C, D);

    //Russell`s approximation
    cout << Russell(S, C, D);

    return 0;
}
