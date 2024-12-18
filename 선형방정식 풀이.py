import re

def parse_equations(equations):
    """선형 방정식을 파싱하여 계수 행렬(A)과 상수 벡터(b)로 변환"""
    variables = []
    parsed_equations = []

    # 방정식 분석
    for eq in equations:
        terms = re.findall(r'([+-]?\s*\d*\.?\d*)([a-zA-Z]+)', eq)
        constant = re.findall(r'[=]\s*([+-]?\s*\d*\.?\d*)', eq)
        eq_dict = {}

        for coef, var in terms:
            if var not in variables:
                variables.append(var)
            coef = coef.replace(' ', '')
            eq_dict[var] = float(coef) if coef not in ('', '+', '-') else float(f"{coef}1")

        # 상수 추출
        eq_dict['constant'] = -float(constant[0]) if constant else 0
        parsed_equations.append(eq_dict)

    # 변수 정렬
    variables = sorted(variables)

    # 계수 행렬(A)와 상수 벡터(b) 생성
    A = []
    b = []

    for eq_dict in parsed_equations:
        row = [eq_dict.get(var, 0) for var in variables]
        A.append(row)
        b.append([-eq_dict['constant']])

    return A, b, variables

def determinant(matrix):
    """행렬식 계산 (재귀적으로 구현)"""
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

    det = 0
    for c in range(n):
        sub_matrix = [row[:c] + row[c+1:] for row in matrix[1:]]
        det += ((-1) ** c) * matrix[0][c] * determinant(sub_matrix)
    return det

def inverse(matrix):
    """행렬의 역행렬 계산"""
    n = len(matrix)
    det = determinant(matrix)
    if det == 0:
        return None

    # 여인자 행렬
    cofactors = []
    for r in range(n):
        cofactor_row = []
        for c in range(n):
            minor = [row[:c] + row[c+1:] for row in (matrix[:r] + matrix[r+1:])]
            cofactor_row.append(((-1) ** (r + c)) * determinant(minor))
        cofactors.append(cofactor_row)

    # 전치 및 스케일링
    cofactors = list(map(list, zip(*cofactors)))  # 전치
    for r in range(n):
        for c in range(n):
            cofactors[r][c] /= det
    return cofactors

def matrix_mult(A, B):
    """두 행렬의 곱 계산"""
    return [[sum(a * b for a, b in zip(row, col))] for row in A for col in zip(*B)]

def gauss_jordan_with_free_vars(matrix, b):
    """가우스-조던 소거법 및 자유변수 도입"""
    n = len(matrix)
    augmented = [matrix[i] + [b[i][0]] for i in range(n)]  # 증강행렬 생성

    # 전진 소거 단계
    for i in range(n):
        if augmented[i][i] == 0:
            for j in range(i + 1, n):
                if augmented[j][i] != 0:
                    augmented[i], augmented[j] = augmented[j], augmented[i]
                    break
        if augmented[i][i] == 0:
            continue

        pivot = augmented[i][i]
        augmented[i] = [x / pivot for x in augmented[i]]

        for j in range(n):
            if i != j:
                factor = augmented[j][i]
                augmented[j] = [augmented[j][k] - factor * augmented[i][k] for k in range(len(augmented[0]))]

    rank_A = sum(1 for row in augmented if any(abs(x) > 1e-9 for x in row[:-1]))
    rank_augmented = sum(1 for row in augmented if any(abs(x) > 1e-9 for x in row))

    if rank_A < rank_augmented:
        return "불능"
    elif rank_A < len(matrix):
        free_vars = []
        solutions = [None] * n
        for i in range(n):
            if all(abs(augmented[i][j]) < 1e-9 for j in range(len(augmented[i]) - 1)):
                free_vars.append(i)
            else:
                solutions[i] = augmented[i][-1]

        general_solution = []
        free_var_idx = 0
        for i in range(n):
            if i in free_vars:
                general_solution.append(f"t{free_var_idx}")
                free_var_idx += 1
            else:
                equation = solutions[i]
                for j in free_vars:
                    equation -= augmented[i][j] * f"t{free_vars.index(j)}"
                general_solution.append(equation)
        return general_solution
    else:
        return [[row[-1]] for row in augmented]

def solve_linear_system():
    """메인 함수"""
    print("선형 연립방정식 풀이기 (예: '2x + 3y = 8', '-x + 4y = 3')")
    n = int(input("방정식의 개수를 입력하세요: "))
    
    print("방정식을 입력하세요:")
    equations = [input(f"{i+1}번째 방정식: ") for i in range(n)]

    A, b, variables = parse_equations(equations)

    print("\n입력된 A 행렬:")
    for row in A:
        print(row)
    print("입력된 b 벡터:")
    for row in b:
        print(row)

    print("\n변수 순서:", variables)

    if len(A) == len(A[0]):
        det = determinant(A)
        if det != 0:
            print("\nA는 가역행렬입니다. 역행렬을 사용하여 해를 구합니다.")
            A_inv = inverse(A)
            if A_inv is None:
                print("역행렬을 계산할 수 없습니다.")
                return
            x = matrix_mult(A_inv, b)
            print("해 x는 다음과 같습니다:")
            for i, value in enumerate(x):
                print(f"{variables[i]} = {value[0]}")
        else:
            print("\nA는 비가역행렬입니다. 가우스-조던 소거법을 사용합니다.")
            solution = gauss_jordan_with_free_vars(A, b)
            if solution == "불능":
                print("해가 존재하지 않습니다. 불능입니다.")
            elif isinstance(solution, list):
                print("해가 무한히 많습니다. 자유변수를 도입한 일반해는 다음과 같습니다:")
                for i, sol in enumerate(solution):
                    print(f"{variables[i]} = {sol}")
            else:
                print("특정 해를 구했습니다:")
                for i, value in enumerate(solution):
                    print(f"{variables[i]} = {value[0]}")
    else:
        print("A가 정사각행렬이 아닙니다. 해를 구할 수 없습니다.")

# 실행
solve_linear_system()
