#Arg_list:F-field,n-number of blocks,r-number of rounds
#Create affine layers(matrix and round constant) for LowMC
def create_affine_layers(F,n,r):
    Encrypt_mat_array=[]
    Keygen_mat_array=[]
    Round_constant_array=[]
    for i in range(r):
        Encrypt_mat=Matrix(F,n,lambda i,j:F.random_element())
        Keygen_mat=Matrix(F,n,lambda i,j:F.random_element())
        while rank(Encrypt_mat)<n:
            Encrypt_mat=Matrix(F,n,lambda i,j:F.random_element())
        while rank(Keygen_mat)<n:
            Keygen_mat=Matrix(F,n,lambda i,j:F.random_element())
        Encrypt_mat_array.append(Encrypt_mat)
        Keygen_mat_array.append(Keygen_mat)
        random_vec=[0 for _ in range(n)]
        for j in range(n):
            random_vec[j]=F.random_element()
        Round_constant=vector(random_vec)
        Round_constant_array.append(Round_constant)
    return Encrypt_mat_array,Keygen_mat_array,Round_constant_array

#Arg_list:mat-matrix for affine layers,constant-round constant for affine layers,v-bit vector
#Affine transformation for LowMC
def affine_layer(mat,constant,v):
    return mat*v+constant

#Arg_list:mat-matrix for affine layers,constant-round constant for affine layers,v-bit vector
#Inverse affine transformation for LowMC
def affine_layer_inverse(mat,constant,v):
    return mat.inverse()*(v+constant)

#Arg_list:v-bit vector
#Sbox for LowMC
def Sbox(v):
    result=[]
    y0=v[0]+v[1]*v[2]
    y1=v[0]+v[1]+v[0]*v[2]
    y2=v[0]+v[1]+v[2]+v[0]*v[1]
    result.append(y0)
    result.append(y1)
    result.append(y2)
    result=vector(result)
    return result

#Arg_list:v-bit vector
#Inverse sbox for LowMC
def Sbox_inverse(v):
    result=[]
    y0=v[0]+v[1]+v[1]*v[2]
    y1=v[1]+v[0]*v[2]
    y2=v[0]+v[1]+v[2]+v[0]*v[1]
    result.append(y0)
    result.append(y1)
    result.append(y2)
    result=vector(result)
    return result

#Arg_list:v-bit vector
#Nolinear transformation for LowMC
def S(v):
    result=[]
    for i in range(n//3):
        tmp=vector(v[i*3:(i+1)*3])
        result+=(Sbox(tmp))
    result=vector(result)
    return result

#Arg_list:v-bit vector
#Inverse nolinear transformation for LowMC
def S_inverse(v):
    result=[]
    for i in range(n//3):
        tmp=vector(v[i*3:(i+1)*3])
        result+=(Sbox_inverse(tmp))
    result=vector(result)
    return result

#Arg_list:v-bit vector,key-secret key for encryption,encrypt_mat_array,keygen_mat_array-matrix array,round_constant_array-round constant array,start_r-start round for encryption,end_r-end round for encryption
#Encrypt a bit vector v to its ciphertext
def LowMC(v,key,Encrypt_mat_array,Keygen_mat_array,Round_constant_array,start_r,end_r):
    v+=key
    for i in range(start_r,end_r):
        v=S(v)
        v=affine_layer(Encrypt_mat_array[i],Round_constant_array[i],v)
        v+=(Keygen_mat_array[i]*key)
    return v

#Arg_list:v-bit vector,key-secret key for encryption,encrypt_mat_array,keygen_mat_array-matrix array,round_constant_array-round constant array,start_r-start round for encryption,end_r-end round for encryption
#Dncrypt a bit vector v to its plaintext
def LowMC_inverse(v,key,Encrypt_mat_array,Keygen_mat_array,Round_constant_array,start_r,end_r):
    for i in range(start_r,end_r,-1):
        v+=(Keygen_mat_array[i]*key)
        v=affine_layer_inverse(Encrypt_mat_array[i],Round_constant_array[i],v)
        v=S_inverse(v)
    return v

#Arg_list:F-field,n-number of blocks,r-number of rounds,key-secret key for encryption,encrypt_mat_array,keygen_mat_array-matrix array,round_constant_array-round constant array
#Generate polynomials
def poly_gen(F,n,r,k,key,encrypt_mat_array,keygen_mat_array,round_constant_array):
    #Key vector for polynomial
    key_v=[F(1) for _ in range(n)]
    #Guess 2 input bits  of Sbox and set 1 input bit as variant
    for i in range(n//3):
    #     key_v[3*i+0]=F.random_element()
    #     key_v[3*i+1]=F.random_element()
        key_v[3*i+2]=k[i]
    key_v=vector(key_v)
    key=vector([F(1) for _ in range(n)])

    #The first plaintext and ciphertext pair
    plain0=vector([F(0) for _ in range(n)])
    for i in range(n):
        plain0[i]=F.random_element()
    cipher0=LowMC(plain0,key,encrypt_mat_array,keygen_mat_array,round_constant_array,0,r)

    #Differences
    deta=vector([F(0) for _ in range(n)])
    deta[n-1]=F(1)
#     for i in range(n//3):
#         deta[3*i+2]=F(1)

    #The second plaintext and ciphertext pair
    plain1=plain0+deta
    cipher1=LowMC(plain1,key,encrypt_mat_array,keygen_mat_array,round_constant_array,0,r)

    #Polynomials from the first plaintext and ciphertext pair
    poly0=[]
    internal_state00=LowMC(plain0,key_v,encrypt_mat_array,keygen_mat_array,round_constant_array,0,ceil(r/2))
    internal_state01=LowMC_inverse(cipher0,key_v,encrypt_mat_array,keygen_mat_array,round_constant_array,r-1,floor(r/2))
    for i in range(n):
        poly0.append(internal_state00[i]-internal_state01[i])
    
    #Polynomials from the second plaintext and ciphertext pair
    poly1=[]
    internal_state10=LowMC(plain1,key_v,encrypt_mat_array,keygen_mat_array,round_constant_array,0,ceil(r/2))
    internal_state11=LowMC_inverse(cipher1,key_v,encrypt_mat_array,keygen_mat_array,round_constant_array,r-1,floor(r/2))
    for i in range(n):
        poly1.append(internal_state10[i]-internal_state11[i])

    #Final polynomials
    poly=[]
    for i in range(n):
        poly.append(poly0[i]+poly1[i])
    return poly

#Test if polynomials are correct or not
def test_poly(poly,key):
    root=[0 for _ in range(n/3)]
    for i in range(len(root)):
        root[i]=key[3*i+2]
    for i in range(len(poly)):
        if poly[i](root)!=0:
            print("False:poly{} is incorrect!".format(i))
            return 0
    print("All polynomials check passed")
    return 1

def print_poly(poly):
    for i in range(len(poly)):
        print("poly{}:={}:".format(i,poly[i]))

if __name__ == "__main__":
    #Arg generation
    F=GF(2)#field
    n=21  #blocksize
    r=7    #number of rounds
    K=PolynomialRing(F,n/3,'k') #Polynomial ring
    k=K.gens()  #Variates
    encrypt_mat_array,keygen_mat_array,round_constant_array=create_affine_layers(F,n,r)
    key=vector([F(1) for _ in range(n)])
    #Polynomials generation
    poly=poly_gen(F,n,r,k,key,encrypt_mat_array,keygen_mat_array,round_constant_array)   
    print("Polynomials have all been generated")
    #Output polynomials to file
    with open(f"./experiments/system_{r}_round.txt", 'w') as f:
        for polynomial in poly:
            f.write(f"{polynomial}\n")
#     print_poly(poly)
    test_poly(poly,key)