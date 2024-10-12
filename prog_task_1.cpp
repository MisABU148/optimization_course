#include <iostream>

using namespace std;

int n = 4, m = 8;
double A[4][8] {};
bool is_appl = 1;

bool Simplex_method() {
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

int main() {
    double c[6];
    double b[3];
    for (int i = 0; i < 3;i++) {
        cin >> c[i];
        A[0][i+1] = -c[i];
    }

    for (int i = 1; i < n;i++) {
        for (int j = 1; j < m-1;j++) {
            cin >> A[i][j];
        }
    }

    for (int i = 0; i < n-1;i++) {
        cin >> b[i];
        A[i+1][m-1] = b[i];
        A[i+1][0] = i+4;
    }

    while (!Simplex_method()) {}

    if (is_appl) {
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
        cout << endl << A[0][m-1];
    } else {
        cout << "The method is not applicable!";
    }

    return 0;
}
