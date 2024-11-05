#include <iostream>
#include <vector>

using namespace std;

int Vogel_iteration(vector<int>& S, vector<vector<int>>& C, vector<int>& D) {
    int s = S.size(), d = D.size();
    int max_dif = 0, min_row = -1, min_col = -1;
    for (int i = 0; i<s; i++) {
        int mn_1 = 1000, mn_2 = 1000;
        for (int j = 0; j < d; j++) {
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

    for (int j = 0; j<d; j++) {
        int mn_1 = 1000, mn_2 = 1000;
        for (int i = 0; i < s; i++) {
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

    int min_weight = 1000;
    if (min_col == -1) {
        for (int i = 0; i<d; i++) {
            if (min_weight > C[min_row][i]) {
                min_col = i;
                min_weight = C[min_row][i];
            }
        }
    } else {
        for (int j = 0; j<s; j++) {
            if (min_weight > C[j][min_col]) {
                min_row = j;
                min_weight = C[j][min_col];
            }
        }
    }

    int way = min(D[min_col], S[min_row]) * C[min_row][min_col];
    if (D[min_col] < S[min_row]) {
        for (auto& row : C) {
            row.erase(row.begin() + min_col);
        }
        S[min_row] -= min(D[min_col], S[min_row]);
        D.erase(D.begin() + min_col);
    } else if (S[min_row] < D[min_col]) {
        D[min_col] -= min(D[min_col], S[min_row]);
        C.erase(C.begin() + min_row);
        S.erase(S.begin() + min_row);
    } else {
        for (auto& row : C) {
            row.erase(row.begin() + min_col);
        }
        D.erase(D.begin() + min_col);
        C.erase(C.begin() + min_row);
        S.erase(S.begin() + min_row);
    }

    return way;
}

int Vogel(vector<int> S, vector<vector<int>> C, vector<int> D) {
    int initial = 0;
    while (S.size()!=0 && D.size()!=0) {
        initial += Vogel_iteration(S, C, D);
    }
    return initial;
}

int main() {
    int s = 3, d = 4;
    int c;
    vector<int> S;
    for (int i = 0; i < s; i++) {
        cin >> c;
        S.push_back(c);
    }
    vector<vector<int>> C;
    for (int i = 0; i < s; i++) {
        C.push_back({});
        for (int j = 0; j < d; j++) {
            cin >> c;
            C[i].push_back(c);
        }
    }
    vector<int> D;
    for (int i = 0; i < d; i++) {
        cin >> c;
        D.push_back(c);
    }
    Vogel(S, C, D);
}
