from numpy import zeros
entradas = open('saved_data.txt', 'r')
x = entradas.readlines()
str_nx = int(float(x[9]))
str_ny = int(float(x[10]))
str_nelx = int(float(x[2]))
str_nely = int(float(x[3]))
n_elem = str_nelx*str_nely
nndes = 2 * str_nx * str_ny
n_gdl = int(float(nndes))
vazia = zeros((int(float(nndes)), int(float(nndes))), dtype=float) #matriz de zeros
entradas.close()
matriz_ele = open('matriz R do elemento.txt', 'r')
a = matriz_ele.readlines()
for l in range(64):
    a[l] = float(a[l])

matriz_ele.close()
kel0 = []
altera = []
i = 0
for l in range(8):
    for c in range(8):
        altera.append(a[i])
        i += 1
    kel0.append(altera[:])
    altera.clear()

a = -2
c = 0
b = -2
x = [0.2, 0.3, 0.5, 0.7]
p = float(3)
for i in range(n_elem):
    a+=2
    b+=2
    for l in range(8):
        if (l == 0):
            vazia[b][a] += kel0[0][c]*(x[i]**p)
            vazia[b][a + 1] += kel0[0][c + 1]*(x[i]**p)
            vazia[b][a + (2 * str_ny)] += kel0[0][c + 2]*(x[i]**p)
            vazia[b][a + (2 * str_ny) + 1] += kel0[0][c + 3]*(x[i]**p)
            vazia[b][a + (2 * str_ny) + 2] += kel0[0][c + 4]*(x[i]**p)
            vazia[b][a + (2 * str_ny) + 3] += kel0[0][c + 5]*(x[i]**p)
            vazia[b][a + 2] += kel0[0][c + 6]*(x[i]**p)
            vazia[b][a + 3] += kel0[0][c + 7]*(x[i]**p)

        elif (l == 1):
            vazia[b + 1][a] += kel0[1][c]*(x[i]**p)
            vazia[b + 1][a + 1] += kel0[1][c + 1]*(x[i]**p)
            vazia[b + 1][a + (2 * str_ny)] += kel0[1][c + 2]*(x[i]**p)
            vazia[b + 1][a + (2 * str_ny) + 1] += kel0[1][c + 3]*(x[i]**p)
            vazia[b + 1][a + (2 * str_ny) + 2] += kel0[1][c + 4]*(x[i]**p)
            vazia[b + 1][a + (2 * str_ny) + 3] += kel0[1][c + 5]*(x[i]**p)
            vazia[b + 1][a + 2] += kel0[1][c + 6]*(x[i]**p)
            vazia[b + 1][a + 3] += kel0[1][c + 7]*(x[i]**p)

        elif (l == 2):
            vazia[b + str_ny * 2][a] += kel0[2][c]*(x[i]**p)
            vazia[b + str_ny * 2][a + 1] += kel0[2][c + 1]*(x[i]**p)
            vazia[b + str_ny * 2][a + (2 * str_ny)] += kel0[2][c + 2]*(x[i]**p)
            vazia[b + str_ny * 2][a + (2 * str_ny) + 1] += kel0[2][c + 3]*(x[i]**p)
            vazia[b + str_ny * 2][a + (2 * str_ny) + 2] += kel0[2][c + 4]*(x[i]**p)
            vazia[b + str_ny * 2][a + (2 * str_ny) + 3] += kel0[2][c + 5]*(x[i]**p)
            vazia[b + str_ny * 2][a + 2] += kel0[2][c + 6]*(x[i]**p)
            vazia[b + str_ny * 2][a + 3] += kel0[2][c + 7]*(x[i]**p)

        elif (l == 3):
            vazia[b + str_ny * 2 + 1][a] += kel0[3][c]*(x[i]**p)
            vazia[b + str_ny * 2 + 1][a + 1] += kel0[3][c + 1]*(x[i]**p)
            vazia[b + str_ny * 2 + 1][a + (2 * str_ny)] += kel0[3][c + 2]*(x[i]**p)
            vazia[b + str_ny * 2 + 1][a + (2 * str_ny) + 1] += kel0[3][c + 3]*(x[i]**p)
            vazia[b + str_ny * 2 + 1][a + (2 * str_ny) + 2] += kel0[3][c + 4]*(x[i]**p)
            vazia[b + str_ny * 2 + 1][a + (2 * str_ny) + 3] += kel0[3][c + 5]*(x[i]**p)
            vazia[b + str_ny * 2 + 1][a + 2] += kel0[3][c + 6]*(x[i]**p)
            vazia[b + str_ny * 2 + 1][a + 3] += kel0[3][c + 7]*(x[i]**p)

        elif (l == 4):
            vazia[b + str_ny * 2 + 2][a] += kel0[4][c]*(x[i]**p)
            vazia[b + str_ny * 2 + 2][a + 1] += kel0[4][c + 1]*(x[i]**p)
            vazia[b + str_ny * 2 + 2][a + (2 * str_ny)] += kel0[4][c + 2]*(x[i]**p)
            vazia[b + str_ny * 2 + 2][a + (2 * str_ny) + 1] += kel0[4][c + 3]*(x[i]**p)
            vazia[b + str_ny * 2 + 2][a + (2 * str_ny) + 2] += kel0[4][c + 4]*(x[i]**p)
            vazia[b + str_ny * 2 + 2][a + (2 * str_ny) + 3] += kel0[4][c + 5]*(x[i]**p)
            vazia[b + str_ny * 2 + 2][a + 2] += kel0[4][c + 6]*(x[i]**p)
            vazia[b + str_ny * 2 + 2][a + 3] += kel0[4][c + 7]*(x[i]**p)

        elif (l == 5):
            vazia[b + str_ny * 2 + 3][a] += kel0[5][c]*(x[i]**p)
            vazia[b + str_ny * 2 + 3][a + 1] += kel0[5][c + 1]*(x[i]**p)
            vazia[b + str_ny * 2 + 3][a + (2 * str_ny)] += kel0[5][c + 2]*(x[i]**p)
            vazia[b + str_ny * 2 + 3][a + (2 * str_ny) + 1] += kel0[5][c + 3]*(x[i]**p)
            vazia[b + str_ny * 2 + 3][a + (2 * str_ny) + 2] += kel0[5][c + 4]*(x[i]**p)
            vazia[b + str_ny * 2 + 3][a + (2 * str_ny) + 3] += kel0[5][c + 5]*(x[i]**p)
            vazia[b + str_ny * 2 + 3][a + 2] += kel0[5][c + 6]*(x[i]**p)
            vazia[b + str_ny * 2 + 3][a + 3] += kel0[5][c + 7]*(x[i]**p)

        elif (l == 6):
            vazia[b + 2][a] += kel0[6][c]*(x[i]**p)
            vazia[b + 2][a + 1] += kel0[6][c + 1]*(x[i]**p)
            vazia[b + 2][a + (2 * str_ny)] += kel0[6][c + 2]*(x[i]**p)
            vazia[b + 2][a + (2 * str_ny) + 1] += kel0[6][c + 3]*(x[i]**p)
            vazia[b + 2][a + (2 * str_ny) + 2] += kel0[6][c + 4]*(x[i]**p)
            vazia[b + 2][a + (2 * str_ny) + 3] += kel0[6][c + 5]*(x[i]**p)
            vazia[b + 2][a + 2] += kel0[6][c + 6]*(x[i]**p)
            vazia[b + 2][a + 3] += kel0[6][c + 7]*(x[i]**p)

        elif (l == 7):
            vazia[b + 3][a] += kel0[7][c]*(x[i]**p)
            vazia[b + 3][a + 1] += kel0[7][c + 1]*(x[i]**p)
            vazia[b + 3][a + (2 * str_ny)] += kel0[7][c + 2]*(x[i]**p)
            vazia[b + 3][a + (2 * str_ny) + 1] += kel0[7][c + 3]*(x[i]**p)
            vazia[b + 3][a + (2 * str_ny) + 2] += kel0[7][c + 4]*(x[i]**p)
            vazia[b + 3][a + (2 * str_ny) + 3] += kel0[7][c + 5]*(x[i]**p)
            vazia[b + 3][a + 2] += kel0[7][c + 6]*(x[i]**p)
            vazia[b + 3][a + 3] += kel0[7][c + 7]*(x[i]**p)


print(f'{vazia}')

file = open('matriz Global.txt', 'w+')

for i in range(n_gdl):
    for j in range(n_gdl):
        file.write(f'{vazia[i][j]}\n')

file.close()


