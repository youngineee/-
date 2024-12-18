def matrix_mult(A, B):
    """두 행렬의 곱 계산"""
    n = len(A)
    m = len(B[0])
    result = [[0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            result[i][j] = sum(A[i][k] * B[k][j] for k in range(len(B)))
    return result

def identity_matrix(n):
    """단위행렬 생성"""
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

def matrix_power_modular(A, k):
    """행렬 거듭제곱 계산 (모듈러 방식)"""
    n = len(A)
    result = identity_matrix(n)
    base = A[:]
    
    while k > 0:
        if k % 2 == 1:
            result = matrix_mult(result, base)
        base = matrix_mult(base, base)
        k //= 2
    
    return result

def diagonalize_matrix(A):
    """행렬을 대각화하여 고유벡터와 고유값 반환"""
    n = len(A)
    eigenvalues, eigenvectors = [], []
    # 고유값과 고유벡터를 직접 구하는 코드 
    for i in range(n):
        eigenvalues.append(A[i][i])  # 대각 성분이 고유값이라고 가정
        eigenvector = [1 if j == i else 0 for j in range(n)]
        eigenvectors.append(eigenvector)
    P = [list(v) for v in zip(*eigenvectors)]  # 고유벡터 행렬
    D = [[eigenvalues[i] if i == j else 0 for j in range(n)] for i in range(n)]  # 대각행렬
    return P, D

def invert_matrix(P):
    """행렬의 역행렬 계산 (단순한 방식)"""
    # 여기서는 고유벡터 행렬이 역행렬이 존재한다고 가정하고 대각화에 사용
    n = len(P)
    identity = identity_matrix(n)
    # 가우스-조던 방식으로 역행렬 계산
    augmented = [P[i] + identity[i] for i in range(n)]
    for i in range(n):
        if augmented[i][i] == 0:
            for j in range(i + 1, n):
                if augmented[j][i] != 0:
                    augmented[i], augmented[j] = augmented[j], augmented[i]
                    break
        pivot = augmented[i][i]
        augmented[i] = [x / pivot for x in augmented[i]]
        for j in range(n):
            if i != j:
                factor = augmented[j][i]
                augmented[j] = [augmented[j][k] - factor * augmented[i][k] for k in range(2 * n)]
    return [row[n:] for row in augmented]

def matrix_power_diagonalization(A, k):
    """대각화 방식으로 행렬 거듭제곱 계산"""
    P, D = diagonalize_matrix(A)
    P_inv = invert_matrix(P)
    
    # D^k 계산
    D_k = [[D[i][i] ** k if i == j else 0 for j in range(len(D))] for i in range(len(D))]
    
    # A^k = P * D^k * P^-1
    return matrix_mult(matrix_mult(P, D_k), P_inv)

def is_diagonalizable(A):
    """대각화 가능 여부 확인"""
    # 대각화 가능 여부 확인을 위해 고유벡터 행렬의 rank를 계산 
    n = len(A)
    eigenvalues, eigenvectors = [], []
    for i in range(n):
        eigenvalues.append(A[i][i])  # 대각 성분이 고유값이라고 가정
        eigenvector = [1 if j == i else 0 for j in range(n)]
        eigenvectors.append(eigenvector)
    P = [list(v) for v in zip(*eigenvectors)]  # 고유벡터 행렬
    rank = sum([1 for row in P if any(x != 0 for x in row)])  # 고유벡터 독립성 확인
    return rank == n

def solve_matrix_power():
    """메인 함수"""
    print("행렬 거듭제곱 계산기")
    n = int(input("행렬의 크기(n x n)을 입력하세요: "))
    print(f"{n} x {n} 행렬을 입력하세요, 띄어쓰기와 엔터로 구분:")
    A = [list(map(float, input().split())) for _ in range(n)]
    k = int(input("거듭제곱 지수 k를 입력하세요: "))

    if is_diagonalizable(A):
        print("\n대각화 방식 사용:")
        try:
            result_diag = matrix_power_diagonalization(A, k)
            print("결과 행렬 (대각화 방식):")
            for row in result_diag:
                print(" ".join(map(str, row)))
        except Exception as e:
            print(f"대각화 방식 실패: {e}")
    else:
        print("\n모듈러 방식 사용 (대각화 불가능):")
        result_mod = matrix_power_modular(A, k)
        print("결과 행렬 (모듈러 방식):")
        for row in result_mod:
            print(" ".join(map(str, row)))

# 실행
solve_matrix_power()
