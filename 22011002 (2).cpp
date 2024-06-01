#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#define SIZE 10
#define TOLERANCE 1e-6
#define MAX_ITERATIONS 100
#define MAXEXPR 100
#define MAX 20
#define MAXSTACK 100

void bisection();
void regulaFalsi();
void newtonRaphson();
void inverseMatrix();
	void print_matrix(double** mat, int n);
	double** allocate_matrix(int n);
	void free_matrix(double** mat, int n);
	double determinant(double** mat, int n);
	void get_cofactor(double** mat, double** temp, int p, int q, int n);
	void adjoint(double** mat, double** adj, int n);
	int invert_matrix(double** mat, double** inv, int n);
void gaussEleminasyon();
	void printMatrix(float **matrix, int size);
	void gaussElimination(float **matrix, int size);
void gaussSeidel();
	void gauss(double A[SIZE][SIZE], double b[SIZE], double x[SIZE], int n);
void numericalDifferantial();
void simpson();
void trapez();
void gregoryNewton();
void makeStringANiceInt();
void inverseMatrix();

typedef enum { NUMBER, VARIABLE, operat, FUNCTION, PARENTHESIS } TokenType;

typedef struct {
    TokenType type;
    union {
        double number1;
        char operat;
        char function[10];
    } value;
} Token;

typedef struct {
    Token data[MAXSTACK];
    int top;
} Stack;

void initStack(Stack *s) {
    s->top = -1;
}

int isEmpty(Stack *s) {
    return s->top == -1;
}

int isFull(Stack *s) {
    return s->top == MAXSTACK - 1;
}

void push(Stack *s, Token val) {
	
    if (!isFull(s)) {
        s->data[++(s->top)] = val;
    }
}

Token pop(Stack *s) {
    if (!isEmpty(s)) {
        return s->data[(s->top)--];
    }
    Token emptyToken = {NUMBER, {0}};
    return emptyToken; // Error case
}

Token peek(Stack *s) {
    if (!isEmpty(s)) {
        return s->data[s->top];
    }
    Token emptyToken = {NUMBER, {0}};
    return emptyToken; // Error case
}

int precedence(char op) {
    switch (op) {
        case '+':
        case '-':
            return 1;
        case '*':
        case '/':
            return 2;
        case '^':
            return 3;
        default:
            return 0;
    }
}

int isLeftAssociative(char op) {
    switch (op) {
        case '+':
        case '-':
        case '*':
        case '/':
            return 1;
        case '^':
            return 0;
        default:
            return 1;
    }
}

int isFunction(char *str) {
    return !strcmp(str, "sin") || !strcmp(str, "cos") || !strcmp(str, "tan") ||
           !strcmp(str, "asin") || !strcmp(str, "acos") || !strcmp(str, "atan") ||
           !strcmp(str, "log") || !strcmp(str, "ln");
}

Token getToken(char *expr, int *index, double xValue) {
    Token token;
    token.type = NUMBER;
    token.value.number1 = 0;

    while (isspace(expr[*index])) (*index)++;
    
    if (isdigit(expr[*index]) || (expr[*index] == '.' && isdigit(expr[*index + 1]))) {
        sscanf(expr + *index, "%lf", &token.value.number1);
        token.type = NUMBER;
        while (isdigit(expr[*index]) || expr[*index] == '.') (*index)++;
    } else if (isalpha(expr[*index])) {
        int i = 0;
        while (isalpha(expr[*index])) {
            token.value.function[i++] = expr[*index];
            (*index)++;
        }
        token.value.function[i] = '\0';
        if (isFunction(token.value.function)) {
            token.type = FUNCTION;
        } else if (strcmp(token.value.function, "x") == 0) {
            token.type = VARIABLE;
            token.value.number1 = xValue;
        } else {
            // Error Case
        }
    } else {
        token.value.operat = expr[*index];
        if (strchr("+-*/^", expr[*index])) {
            token.type = operat;
        } else if (strchr("()", expr[*index])) {
            token.type = PARENTHESIS;
        } else {
            // Error Case
        }
        (*index)++;
    }
    
    return token;
}

void shunting(char *expr, Token outputQueue[], int *outputSize, double xValue) {
    Stack operatStack;
    initStack(&operatStack);

    int index = 0, outIndex = 0;
    while (expr[index] != '\0') {
        Token token = getToken(expr, &index, xValue);
        switch (token.type) {
            case NUMBER:
            case VARIABLE:
                outputQueue[outIndex++] = token;
                break;
            case FUNCTION:
                push(&operatStack, token);
                break;
            case operat:
                while (!isEmpty(&operatStack) && operatStack.data[operatStack.top].type != PARENTHESIS &&
                       (precedence(operatStack.data[operatStack.top].value.operat) > precedence(token.value.operat) ||
                       (precedence(operatStack.data[operatStack.top].value.operat) == precedence(token.value.operat) && isLeftAssociative(token.value.operat)))) {
                    outputQueue[outIndex++] = pop(&operatStack);
                }
                push(&operatStack, token);
                break;
            case PARENTHESIS:
                if (token.value.operat == '(') {
                    push(&operatStack, token);
                } else {
                    while (!isEmpty(&operatStack) && operatStack.data[operatStack.top].value.operat != '(') {
                        outputQueue[outIndex++] = pop(&operatStack);
                    }
                    pop(&operatStack); // Pop '('
                    if (!isEmpty(&operatStack) && operatStack.data[operatStack.top].type == FUNCTION) {
                        outputQueue[outIndex++] = pop(&operatStack);
                    }
                }
                break;
        }
    }
    while (!isEmpty(&operatStack)) {
        outputQueue[outIndex++] = pop(&operatStack);
    }
    *outputSize = outIndex;
}

double applyFunction(char *func, double value) {
    if (strcmp(func, "sin") == 0) return sin(value);
    if (strcmp(func, "cos") == 0) return cos(value);
    if (strcmp(func, "tan") == 0) return tan(value);
    if (strcmp(func, "asin") == 0) return asin(value);
    if (strcmp(func, "acos") == 0) return acos(value);
    if (strcmp(func, "atan") == 0) return atan(value);
    if (strcmp(func, "log") == 0) return log10(value);
    if (strcmp(func, "ln") == 0) return log(value);
    return 0; // Error case
}

double applyoperat(char op, double a, double b) {
    switch (op) {
        case '+': return a + b;
        case '-': return a - b;
        case '*': return a * b;
        case '/': return a / b;
        case '^': return pow(a, b);
        default: return 0; // Error case
    }
}

double evaluatePostfix(Token *postfix, int size, double xValue) {
    Stack stack;
    initStack(&stack);
    for (int i = 0; i < size; i++) {
        Token token = postfix[i];
        if (token.type == NUMBER || token.type == VARIABLE) {
            if (token.type == VARIABLE) {
                token.value.number1 = xValue;
            }
            push(&stack, token);
        } else if (token.type == operat) {
            double b = pop(&stack).value.number1;
            double a = pop(&stack).value.number1;
            double result = applyoperat(token.value.operat, a, b);
            Token resultToken = {NUMBER, {result}};
            push(&stack, resultToken);
        } else if (token.type == FUNCTION) {
            double a = pop(&stack).value.number1;
            double result = applyFunction(token.value.function, a);
            Token resultToken = {NUMBER, {result}};
            push(&stack, resultToken);
        }
    }
    return pop(&stack).value.number1;
}

double f(double x, Token *postfix, int size) {
    return evaluatePostfix(postfix, size, x);
}

int main() {
    int a;
    printf("\n1 bisection yontemi\n2 regula-falsi yontemi\n3 newton-rapshon yontemi\n4 NxN'lik bir matrisin tersi\n5 gauss eleminasyon\n6 gauss seidal yontemleri\n7 sayisal turev\n8 simpson yontemi\n9 trapez yontemi\n10 degisken donusumsuz gregory enterpolasyonu");
    printf("\nYapmak istediginiz yontemi seciniz: ");
    scanf("%d", &a);
    getchar(); // scanf sonrasi giris tamponunu temizlemek icin

    switch(a) { //hangi algoritmanin caliscacaginin secildigi yer
        case 1:
            printf("Simdi bisection\n");//+
            bisection(); 
            break;
        case 2:
            printf("Simdi regula-falsi\n");//+
			regulaFalsi();
            break;
        case 3:
            printf("Simdi newton-rapshon\n");//+
            newtonRaphson();
            break;
        case 4:
            printf("NxN'lik bir matrisin tersi\n");//+
            inverseMatrix();
            break;
        case 5:
            printf("Gauss eleminasyon\n");//+
            gaussEleminasyon();
            break;
        case 6:
            printf("Gauss seidal\n");//+
            gaussSeidel();
            break;
        case 7:
            printf("sayisal turev\n");//+
            numericalDifferantial();
            break;
        case 8:
            printf("Simpson\n");
            simpson();
            break;
        case 9:
            printf("Trapez\n");//+
            trapez();
            break;
        case 10:
            printf("Degisken donusumsuz gregory newton enterpolasyonu\n");//+
            gregoryNewton();//+
            break;
        default:
            printf("Bir seyler ters gitti\n");
    }

    return 0;
}
void bisec(double a, double b, Token *postfix, int size, double tol) {
    double fa = f(a, postfix, size);
    double fb = f(b, postfix, size);

    if (fa * fb >= 0) {
        printf("Baslangic araliginda kok yok.\n");
        return;
    }

    double c, fc;
    while ((b - a) / 2.0 > tol) {
        c = (a + b) / 2.0;
        fc = f(c, postfix, size);
        if (fc == 0.0) {
            break;
        } else if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    printf("Kok yaklasik olarak: %lf\n", (a + b) / 2.0);
}void bisection() {
    char expr[MAXEXPR];
    double a, b, tol;
    printf("Matematiksel ifadeyi girin: ");
    fgets(expr, MAXEXPR, stdin);
    expr[strcspn(expr, "\n")] = 0; 

    printf("Baslangic araligini girin (a ve b): ");
    scanf("%lf %lf", &a, &b);

    printf("Tolerans degerini girin: ");
    scanf("%lf", &tol);
    getchar(); // scanf sonrasi giris tamponunu temizlemek icin gerekli

    Token outputQueue[MAXEXPR];
    int outputSize = 0;

    shunting(expr, outputQueue, &outputSize, 0);
    bisec(a, b, outputQueue, outputSize, tol);
}

void regula(double a, double b, Token *postfix, int size, double tol, int maxIter) {
    double fa = f(a, postfix, size);
    double fb = f(b, postfix, size);

    if (fa * fb >= 0) {
        printf("Baslangýc aralýgýnda kok yok.\n");
        return;
    }

    double c, fc;
    for (int iter = 0; iter < maxIter; iter++) {
        c = (a * fb - b * fa) / (fb - fa); // Regula Falsi formülü
        fc = f(c, postfix, size);

        if (fabs(fc) < tol) {
            printf("Kok bulundu: %lf\n", c);
            return;
        }

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    printf("Maksimum iterasyona ulasýldý. Bulunan en iyi tahmin: %lf\n", c);
}void regulaFalsi(){
    char expr[MAXEXPR];
    Token postfix[MAXEXPR];
    int postfixSize;

    printf("Fonksiyonu girin : ");
    fgets(expr, MAXEXPR, stdin);
    expr[strcspn(expr, "\n")] = 0; // fgets kullanýnca gelen \n karakterini silmek için

    shunting(expr, postfix, &postfixSize, 0);

    double a, b;
    double tol;
    int maxIter;

    printf("Baslangic araligi (a ve b) girin: ");
    scanf("%lf %lf", &a, &b);

    printf("Hata toleransini girin: ");
    scanf("%lf", &tol);

    printf("Maksimum iterasyon sayisini girin: ");
    scanf("%d", &maxIter);

    regula(a, b, postfix, postfixSize, tol, maxIter);

}
void newton(double initialGuess, Token *postfix, int size, Token *postfixDeriv, int sizeDeriv, double tol, int maxIter) {
    double x0 = initialGuess;
    double x1, f0, f0_deriv;

    for (int iter = 0; iter < maxIter; iter++) {
        f0 = f(x0, postfix, size);
        f0_deriv = f(x0, postfixDeriv, sizeDeriv);

        if (fabs(f0_deriv) < 1e-12) {
            printf("Turev çok kucuk, Newton-Raphson yontemi durduruldu.\n");
            return;
        }

        x1 = x0 - f0 / f0_deriv;

        printf("Iterasyon %d: x = %lf, f(x) = %lf\n", iter + 1, x1, f0);

        if (fabs(x1 - x0) < tol) {
            printf("Kok bulundu: %lf\n", x1);
            return;
        }

        x0 = x1;
    }

    printf("Maksimum iterasyon sayisina ulasildi. Sonuc: %lf\n", x1);
}void newtonRaphson(){
	char expr[MAXEXPR];
    char exprDeriv[MAXEXPR];
    Token postfix[MAXEXPR];
    Token postfixDeriv[MAXEXPR];
    int postfixSize, postfixSizeDeriv;

    printf("Fonksiyonu giriniz: ");
    fgets(expr, MAXEXPR, stdin);
    expr[strcspn(expr, "\n")] = 0; // fgets kullaninca gelen \n karakterini silmek için

    printf("Fonksiyonun turevini giriniz: ");
    fgets(exprDeriv, MAXEXPR, stdin);
    exprDeriv[strcspn(exprDeriv, "\n")] = 0; // fgets kullaninca gelen \n karakterini silmek için

    shunting(expr, postfix, &postfixSize, 0);
    shunting(exprDeriv, postfixDeriv, &postfixSizeDeriv, 0);

    double initialGuess;
    double tol;
    int maxIter;

    printf("Baslangic tahminini girin: ");
    scanf("%lf", &initialGuess);

    printf("Hata toleransini girin: ");
    scanf("%lf", &tol);

    printf("Maksimum iterasyon sayisini girin: ");
    scanf("%d", &maxIter);

    newton(initialGuess, postfix, postfixSize, postfixDeriv, postfixSizeDeriv, tol, maxIter);
}
void inverseMatrix(){

    int n, i, j; // Matris boyutunu buradan ayarlayabilirsiniz

    printf("Matrisinizin kac satirdan olusacagini giriniz: ");
    scanf("%d", &n);

    double** mat = allocate_matrix(n);
    double** inv = allocate_matrix(n);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("[%d][%d] matrisinin degerini giriniz: ", i, j);
            scanf("%lf", &mat[i][j]);
        }
    }

    if (invert_matrix(mat, inv, n)) {
        printf("Matrisin tersi:\n");
        print_matrix(inv, n);
    } else {
        printf("Matrisin tersi alinamaz.\n");
    }

    free_matrix(mat, n);
    free_matrix(inv, n);
}void print_matrix(double** mat, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", mat[i][j]);
        }
        printf("\n");
    }
}double** allocate_matrix(int n) {
    double** mat = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        mat[i] = (double*)malloc(n * sizeof(double));
    }
    return mat;
}void free_matrix(double** mat, int n) {
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}double determinant(double** mat, int n) {
    if (n == 1)
        return mat[0][0];

    double det = 0;
    double** temp = allocate_matrix(n);
    int sign = 1;

    for (int f = 0; f < n; f++) {
        get_cofactor(mat, temp, 0, f, n);
        det += sign * mat[0][f] * determinant(temp, n - 1);
        sign = -sign;
    }

    free_matrix(temp, n);
    return det;
}void get_cofactor(double** mat, double** temp, int p, int q, int n) {
    int i = 0, j = 0;

    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                temp[i][j++] = mat[row][col];

                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}void adjoint(double** mat, double** adj, int n) {
    if (n == 1) {
        adj[0][0] = 1;
        return;
    }

    int sign = 1;
    double** temp = allocate_matrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            get_cofactor(mat, temp, i, j, n);

            sign = ((i + j) % 2 == 0) ? 1 : -1;

            adj[j][i] = (sign) * (determinant(temp, n - 1));
        }
    }

    free_matrix(temp, n);
}int invert_matrix(double** mat, double** inv, int n) {
    double det = determinant(mat, n);
    if (det == 0) {
        return 0;
    }

    double** adj = allocate_matrix(n);
    adjoint(mat, adj, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i][j] = adj[i][j] / det;
        }
    }

    free_matrix(adj, n);
    return 1;
}
void gaussEleminasyon(){
    int size;
    printf("Matris boyutunu girin: ");
    scanf("%d", &size);

    // Dinamik olarak matris olusturma
    float **matrix = (float **)malloc(size * sizeof(float *));
    for (int i = 0; i < size; i++) {
        matrix[i] = (float *)malloc((size + 1) * sizeof(float));
    }

    printf("Matris elemanlarini girin:\n");
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size + 1; j++) {
            scanf("%f", &matrix[i][j]);
        }
    }

    printf("Baslangic matrisi:\n");
    printMatrix(matrix, size);

    gaussElimination(matrix, size);

    printf("ust ucgen formundaki matris:\n");
    printMatrix(matrix, size);

    // Matris bellegini serbest birakma
    for (int i = 0; i < size; i++) {
        free(matrix[i]);
    }
    free(matrix);
}void printMatrix(float **matrix, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size + 1; j++) {
            printf("%8.3f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}void gaussElimination(float **matrix, int size) {
    for (int i = 0; i < size; i++) {
        // Pivoting
        for (int k = i + 1; k < size; k++) {
            if (matrix[i][i] == 0) {
                printf("Matris cozulemez, cunku pivot sifir.\n");
                return;
            }
            float ratio = matrix[k][i] / matrix[i][i];
            for (int j = i; j < size + 1; j++) {
                matrix[k][j] -= ratio * matrix[i][j];
            }
        }
    }
}
void gaussSeidel(){
    int n;

    printf("Denklem sayisini giriniz: ");
    scanf("%d", &n);

    double A[SIZE][SIZE], b[SIZE], x[SIZE];

    printf("Matrisin degerlerini giriniz A:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }

    printf("Matrisin ikincil degerlerini giriniz b:\n");
    for (int i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }

    gauss(A, b, x, n);

    printf("Denklemlerin sonuclarýný giriniz x:\n");
    for (int i = 0; i < n; i++) {
        printf("%lf\n", x[i]);
    }
    
    

}void gauss(double A[SIZE][SIZE], double b[SIZE], double x[SIZE], int n) {
    double x_old[SIZE];
    double sum;
    int iterations = 0;

    for (int i = 0; i < n; i++) {
        x[i] = 0.0;
    }

    while (iterations < MAX_ITERATIONS) {
        for (int i = 0; i < n; i++) {
            x_old[i] = x[i];
        }

        for (int i = 0; i < n; i++) {
            sum = b[i];
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    sum -= A[i][j] * x[j];
                }
            }
            x[i] = sum / A[i][i];
        }

        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = fabs(x[i] - x_old[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        if (max_diff < TOLERANCE) {
            break;
        }

        iterations++;
    }

    if (iterations == MAX_ITERATIONS) {
        printf("Maximum iterasyona ulaþýldý.\n");
    } else {
        printf("iterasyonlarla %d ulasýldý.\n", iterations);
    }
}
double forwardDifference(double x, double h, Token *postfix, int size) {
    return (f(x + h, postfix, size) - f(x, postfix, size)) / h;
}

double backwardDifference(double x, double h, Token *postfix, int size) {
    return (f(x, postfix, size) - f(x - h, postfix, size)) / h;
}

double centralDifference(double x, double h, Token *postfix, int size) {
    return (f(x + h, postfix, size) - f(x - h, postfix, size)) / (2 * h);
}
void numericalDifferantial(){
    char expr[MAXEXPR];
    double x, h;
    int method;

    printf("Fonksiyon ifadesini girin (örnek: 3*x^2 + 2*x + 1): ");
    fgets(expr, MAXEXPR, stdin);
    expr[strcspn(expr, "\n")] = '\0'; 

    Token postfix[MAXEXPR];
    int postfixSize;
    shunting(expr, postfix, &postfixSize, 0);

    printf("Turev alinacak noktayi girin (x): ");
    scanf("%lf", &x);

    printf("Adým buyuklugunu girin (h): ");
    scanf("%lf", &h);

    printf("Yontemi secin (1 - Ileri Fark, 2 - Geri Fark, 3 - Merkez Fark): ");
    scanf("%d", &method);

    double result;
    switch (method) {
        case 1:
            result = forwardDifference(x, h, postfix, postfixSize);
            printf("Ileri Fark Yontemi Sonuc: %lf\n", result);
            break;
        case 2:
            result = backwardDifference(x, h, postfix, postfixSize);
            printf("Geri Fark Yontemi Sonuc: %lf\n", result);
            break;
        case 3:
            result = centralDifference(x, h, postfix, postfixSize);
            printf("Merkez Fark Yontemi Sonuc: %lf\n", result);
            break;
        default:
            printf("Gecersiz yontem secimi.\n");
            break;
    }
}

void simpson(){
}
double trapezon(double a, double b, int n, Token *postfix, int size) {
    double h = (b - a) / n;
    double sum = f(a, postfix, size) + f(b, postfix, size);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += 2 * f(x, postfix, size);
    }

    return (h / 2) * sum;
}void trapez(){
    char expr[MAXEXPR];
    Token postfix[MAXEXPR];
    int postfixSize;

    printf("Fonksiyonu giriniz: ");
    fgets(expr, MAXEXPR, stdin);
    expr[strcspn(expr, "\n")] = 0; // fgets kullanýnca gelen \n karakterini silmek için

    shunting(expr, postfix, &postfixSize, 0);

    double a, b;
    int n;

    printf("Baslangýc araligi giriniz:(a): ");
    scanf("%lf", &a);

    printf("Bitis araligini giriniz (b): ");
    scanf("%lf", &b);

    printf("Trapez sayiyisini giriniz (n): ");
    scanf("%d", &n);

    double result = trapezon(a, b, n, postfix, postfixSize);
    printf("Sonuç: %lf\n", result);
}

void gregoryNewton(){
    double x[MAX], y[MAX][MAX];
    double value, sum, u;
    int i, j, n;
    printf("veri noktalarinin sayisini giriniz:");
    scanf("%d", &n);
	printf("\n");
    printf("Veriyi giriniz:\n");
    for (i = 0; i < n; i++) {
        printf("x[%d] = ", i);
        scanf("%lf", &x[i]);
    }
    for (i = 0; i < n; i++) {
        printf("y[%d] = ", i);
        scanf("%lf", &y[i][0]);
    }

    for (i = 1; i < n; i++) {
        for (j = 0; j < n - i; j++) {
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }
    }

    printf("\nx'in degerini giriniz':\n");
    scanf("%lf", &value);

    sum = y[0][0];
    u = (value - x[0]) / (x[1] - x[0]);
    for (i = 1; i < n; i++) {
        sum += (u * y[0][i]) / i;
        u *= (value - x[i]) / (x[i + 1] - x[i]);
    }

    printf("\nEnterpolasyon degeri x = %lf is: %lf\n", value, sum);
    printf("\n\n");
}

