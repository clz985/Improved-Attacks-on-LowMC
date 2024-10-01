#LowMC Arg
F=GF(2)#Field
n=129 #Blocksize
r=4   #Number of rounds   
#Look_up table for two equations of active sboxs
#Input differences                        0                                                       1                                                       2                                                       3                                                       4                                                       5                                                       6                                                         7
equ_t=[[([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0]))],
       [([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,0]),vector(F,[0,1,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,0]),vector(F,[0,1,1])],vector(F,[1,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,1,0]),vector(F,[1,0,1])],vector(F,[1,1])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,1,0]),vector(F,[0,1,1])],vector(F,[1,1]))],
       [([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,0]),vector(F,[0,0,1])],vector(F,[1,1])),([vector(F,[1,0,0]),vector(F,[0,1,1])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,1]),vector(F,[1,1,0])],vector(F,[0,1])),([vector(F,[1,0,1]),vector(F,[0,1,1])],vector(F,[1,1]))],
       [([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,0]),vector(F,[0,1,0])],vector(F,[1,0])),([vector(F,[1,0,0]),vector(F,[0,0,1])],vector(F,[0,1])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,1]),vector(F,[0,1,0])],vector(F,[0,1])),([vector(F,[0,0,1]),vector(F,[1,1,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0]))],
       [([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,1,0]),vector(F,[0,0,1])],vector(F,[1,0])),([vector(F,[0,1,0]),vector(F,[1,0,1])],vector(F,[0,1])),([vector(F,[0,0,1]),vector(F,[1,1,0])],vector(F,[1,1])),([vector(F,[1,1,0]),vector(F,[1,0,1])],vector(F,[0,0]))],
       [([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,0]),vector(F,[0,1,0])],vector(F,[0,1])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,0]),vector(F,[0,1,1])],vector(F,[1,1])),([vector(F,[0,1,0]),vector(F,[0,0,1])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,1,0]),vector(F,[0,0,1])],vector(F,[0,1])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0]))],
       [([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,0]),vector(F,[0,0,1])],vector(F,[1,0])),([vector(F,[1,0,0]),vector(F,[0,1,1])],vector(F,[0,1])),([vector(F,[0,1,0]),vector(F,[0,0,1])],vector(F,[1,1])),([vector(F,[0,1,0]),vector(F,[1,0,1])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0]))],
       [([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,0,0]),vector(F,[0,1,0])],vector(F,[1,1])),([vector(F,[1,0,0]),vector(F,[0,0,1])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,1,0]),vector(F,[0,0,1])],vector(F,[0,1])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0])),([vector(F,[1,1,0]),vector(F,[1,0,1])],vector(F,[1,1]))]]

#Look_up table for input states of active sboxs
#Input differences                             0                                                                              1                                                                            2                                                                        3                                                                         4                                                                          5                                                                          6                                                                        7
state_in_t=[[([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0]))],
            [([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,1])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,1,0])],vector(F,[1,0,1])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[1,0,0])],vector(F,[0,1,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,1,0])],vector(F,[1,1,0]))],
            [([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,1,0]),vector(F,[0,0,0])],vector(F,[1,1,0])),([vector(F,[0,0,0]),vector(F,[0,1,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,1,0]),vector(F,[0,0,0])],vector(F,[1,0,1])),([vector(F,[0,0,0]),vector(F,[0,1,0]),vector(F,[0,0,0])],vector(F,[0,0,1]))],
            [([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,1]),vector(F,[0,0,1])],vector(F,[1,0,1])),([vector(F,[0,0,0]),vector(F,[0,1,0]),vector(F,[0,1,0])],vector(F,[0,0,1])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[1,0,0]),vector(F,[1,0,0])],vector(F,[1,1,1])),([vector(F,[0,0,0]),vector(F,[1,0,0]),vector(F,[1,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0]))],
            [([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[1,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[1,1,1])),([vector(F,[1,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,1])),([vector(F,[1,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,1,0])),([vector(F,[1,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0]))],
            [([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,1]),vector(F,[0,0,0]),vector(F,[0,0,1])],vector(F,[1,1,1])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,1,0]),vector(F,[0,0,0]),vector(F,[0,1,0])],vector(F,[1,1,0])),([vector(F,[1,0,0]),vector(F,[0,0,0]),vector(F,[1,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[1,0,0]),vector(F,[0,0,0]),vector(F,[1,0,0])],vector(F,[0,0,1])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0]))],
            [([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,1,0]),vector(F,[0,1,0]),vector(F,[0,0,0])],vector(F,[1,0,1])),([vector(F,[0,1,0]),vector(F,[0,1,0]),vector(F,[0,0,0])],vector(F,[0,0,1])),([vector(F,[1,0,0]),vector(F,[1,0,0]),vector(F,[0,0,0])],vector(F,[0,1,0])),([vector(F,[1,0,0]),vector(F,[1,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0]))],
            [([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,1]),vector(F,[0,0,1]),vector(F,[0,0,1])],vector(F,[0,1,1])),([vector(F,[0,1,0]),vector(F,[0,1,0]),vector(F,[0,1,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[1,0,0]),vector(F,[1,0,0]),vector(F,[1,0,0])],vector(F,[0,0,1])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[0,0,0]),vector(F,[0,0,0]),vector(F,[0,0,0])],vector(F,[0,0,0])),([vector(F,[1,0,0]),vector(F,[1,0,0]),vector(F,[1,0,0])],vector(F,[0,1,0]))]]

#Arg_list:v-bit vector
#Transform a bit vector v to its index number,such as [0,1,0] to 2
def vector_to_index(v):
    return Integer(v[0])*2^2+Integer(v[1])*2+Integer(v[2])

#Arg_list:x-index number of bit vector,n-length of bit vector
#Transform a number to its bit vector,such as 2 to [0,1,0]
def index_to_vector(x,n):
    result=[0 for _ in range(n)]
    p=2
    for i in range(n):
        remainder=x%p
        result[n-i-1]=remainder
        x=x//p
    return result

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

#Arg_list:F-field,path-file path
#Load affine layers and plaintext from file
def load_affine_layers(F,path):
    encrypt_mat_array=[[[0 for _ in range(n)] for _ in range(n)] for _ in range(r)]
    keygen_mat_array=[[[0 for _ in range(n)] for _ in range(n)] for _ in range(r)]
    round_constant_array=[[0 for _ in range(n)] for _ in range(r)]
    plain0=[0 for _ in range(n)]
    with open(path, 'r') as f:
        content=f.read()
        char_cnt=-1
        for i in range(r):
            char_cnt+=len("encrypt_mati\n")#mat_name
            for j in range(n):
                char_cnt+=len("\n[") #'\n['
                for k in range(n):
                    encrypt_mat_array[i][j][k]=F(content[char_cnt])
                    char_cnt+=2
            encrypt_mat_array[i]=Matrix(encrypt_mat_array[i])
        for i in range(r):
            char_cnt+=len("keygen_mati\n")#mat_name
            for j in range(n):
                char_cnt+=len("\n[") #'\n['
                for k in range(n):
                    keygen_mat_array[i][j][k]=F(content[char_cnt])
                    char_cnt+=2
            keygen_mat_array[i]=Matrix(keygen_mat_array[i])
        char_cnt+=len("\n")
        for i in range(r):
            char_cnt+=len("round_constanti\n(")#round_constant_name
            for j in range(n):
                round_constant_array[i][j]=F(content[char_cnt])
                char_cnt+=3
            round_constant_array[i]=vector(round_constant_array[i])
        char_cnt+=len("plaini\n(")
        for j in range(n):
            plain0[j]=F(content[char_cnt])
            char_cnt+=3
        plain0=vector(plain0)
    return encrypt_mat_array,keygen_mat_array,round_constant_array,plain0

#Arg_list:encrypt_mat_array,keygen_mat_array-matrix array,tound_constant_array-round constant array,plain0-plaintext,path-file path
#Save affine layers and plaintext which are valid for our attack to a file
def save_affine_layers(encrypt_mat_array,keygen_mat_array,round_constant_array,plain0,path):
    with open(path, 'w') as f:
        for i in range(len(encrypt_mat_array)):
            f.write(f"encrypt_mat{i}\n")
            f.write(f"{encrypt_mat_array[i]}\n")
        for i in range(len(keygen_mat_array)):
            f.write(f"keygen_mat{i}\n")
            f.write(f"{keygen_mat_array[i]}\n")
        for i in range(len(round_constant_array)):
            f.write(f"round_constant{i}\n")
            f.write(f"{round_constant_array[i]}\n")
        f.write("plain0\n")
        f.write(f"{plain0}\n")

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

#Arg_list:v-bit vector,key-secret key for encryption,Encrypt_mat_array,Keygen_mat_array-matrix array,Round_constant_array-round constant array,start_r-start round for encryption,end_r-end round for encryption
#Encrypt a bit vector v to its ciphertext
def LowMC(v,key,Encrypt_mat_array,Keygen_mat_array,Round_constant_array,start_r,end_r):
    v+=key
    for i in range(start_r,end_r):
        v=S(v)
        v=affine_layer(Encrypt_mat_array[i],Round_constant_array[i],v)
        v+=(Keygen_mat_array[i]*key)
    return v

#Arg_list:v-bit vector,key-secret key for encryption,Encrypt_mat_array,Keygen_mat_array-matrix array,Round_constant_array-round constant array,start_r-start round for encryption,end_r-end round for encryption
#Decrypt a bit vector v to its plaintext
def LowMC_inverse(v,key,Encrypt_mat_array,Keygen_mat_array,Round_constant_array,start_r,end_r):
    for i in range(start_r,end_r,-1):
        v+=(Keygen_mat_array[i]*key)
        v=affine_layer_inverse(Encrypt_mat_array[i],Round_constant_array[i],v)
        v=S_inverse(v)
    return v

#Arg_list:v1,v2-two input bit vectors,key-secret key for encryption,Encrypt_mat_array-matrix array,Round_constant_array-round constant array,start_r-start round for encryption,end_r-end round for encryption
#Get input and output differences of every sbox for every round 
def LowMC_diff_path(v1,v2,key,Encrypt_mat_array,Keygen_mat_array,Round_constant_array,start_r,end_r):
    Sbox_in_diff=[]
    Sbox_out_diff=[]
    v1+=key
    v2+=key
    active_cnt_list=[]
    for i in range(start_r,end_r):
        Sbox_in_diff.append(v1+v2)
        v1=S(v1)
        v2=S(v2)
        active_cnt=0
        for j in range(n//3):
            if Sbox_in_diff[i][j*3+0]!=0 or Sbox_in_diff[i][j*3+1]!=0 or Sbox_in_diff[i][j*3+2]!=0:
                active_cnt+=1
        active_cnt_list.append(active_cnt)
        Sbox_out_diff.append(v1+v2)
        v1=affine_layer(Encrypt_mat_array[i],Round_constant_array[i],v1)
        v2=affine_layer(Encrypt_mat_array[i],Round_constant_array[i],v2)
        v1+=(Keygen_mat_array[i]*key)
        v2+=(Keygen_mat_array[i]*key)
    return Sbox_in_diff,Sbox_out_diff,active_cnt_list

#Arg_list:v-input bit vectors,key-secret key for encryption,Encrypt_mat_array-matrix array,Round_constant_array-round constant array,start_r-start round for encryption,end_r-end round for encryption
#Get input and output states of every sbox for every round 
def LowMC_state_path(v,key,Encrypt_mat_array,Keygen_mat_array,Round_constant_array,start_r,end_r):
    Sbox_in=[]
    Sbox_out=[]
    v+=key
    for i in range(start_r,end_r):
        Sbox_in.append(v)
        v=S(v)
        Sbox_out.append(v)
        v=affine_layer(Encrypt_mat_array[i],Round_constant_array[i],v)
        v+=(Keygen_mat_array[i]*key)
    return Sbox_in,Sbox_out

#Arg_list:poly-polynomials,root-roots of polynomials
#Test if polynomials are correct or not
def test_poly(poly,root):
    for i in range(len(poly)):
        if poly[i](root)!=0:
            print("False:poly{} is incorrect!".format(i))
            return 0
    print("All polynomials check passed")
    return 1

#Arg_list:result-answer from gaussian elimination,root-roots of equations
#Test if results are correct or not
def test_result(result,root):
    for i in range(len(result)):
        if result[i]==root[i]:
            continue
        if result[i](root)==root[i]:
            continue
        print("False:result{} is incorrect!".format(i))
        return 0
    print("All results check passed")
    return 1

#Arg_list:F-field,x-variants list,var_cnt-number of variables,poly-polynomials
#Get coefficient matrix of polynomials
def poly_to_affine_mat(F,x,var_cnt,poly):
    result_mat=Matrix(F,[[0 for _ in range(var_cnt)] for _ in range(len(poly))])
    result_constant=vector(F,[0 for _ in range(len(poly))])
    for i in range(len(poly)):
        if poly[i]==0:
            continue
        var_index=0
        poly_list=list(poly[i])
        poly_list_len=len(poly_list)
        if poly_list[poly_list_len-1][1]==1:
            poly_list_len-=1
            result_constant[i]=F(1)
        for j in range(poly_list_len):
            while x[var_index]!=poly_list[j][1] and var_index<var_cnt:
                var_index+=1
            result_mat[i,var_index]=F(1)
    return result_mat,result_constant

#Arg_list:monomial-a monomial in polynomial,var_list-variables list
#Get the index of a monomial under a specific term order
def get_var_index(monomial,var_list):
    var_vec=[1 for _ in range(36)]
    result=[]
    if monomial==1:
        result.append(-1)
        return result
    for i in var_list:
        tmp_var_vec=var_vec.copy()
        tmp_var_vec[i]=0
        if monomial(tmp_var_vec)==0:
            result.append(i)
    return result

#Arg_list:F-field,A-matrix,b-vector,col-number of columns
#Gaussian elimination
def gaussian_elimination(F,A,b,col):
    nc=col
    nr=A.nrows()
    non_free_list=[]
    free_list=[]
    var_dic={}
    for i in range(nc):
        var_dic[i]=0
    for j in range(nc):
        max_id=-1
        for i in range(j,nr):
            if A[i][j]==1:
                max_id=i
                A[[j,i]]=A[[i,j]]
                tmp=b[i]
                b[i]=b[j]
                b[j]=tmp
                break
        if max_id<0:
            continue
        for i in range(j+1,nr):
            if A[i][j]==1:
                A[i]+=A[j]
                b[i]+=b[j]
    for i in range(nr-1,-1,-1):
        max_id=-1
        for j in range(nc):
            if A[i][j]==1:
                max_id=j
                var_dic[j]=1
                if j not in non_free_list:
                    non_free_list.append(j)
                break
        if max_id<0:
            continue
        for j in range(i-1,-1,-1):
            if A[j][max_id]==1:
                A[j]+=A[i]
                b[j]+=b[i]
    bas_mat=Matrix([[x[0] for _ in range(nc+1)] for _ in range(nc)])
    bas_mat-=bas_mat
    for i in range(nc):
        for j in range(nc):
            if var_dic[j]==0:
                if i<nr:
                    bas_mat[i,j]=A[i,j]
                else:
                    bas_mat[i,j]=F(0)
                if j not in free_list:
                    free_list.append(j)
        if i<nr:
            bas_mat[i,nc]=b[i]
        else:
            bas_mat[i,nc]=F(0)
    row_index=0
    non_free_list.reverse()
    for i in non_free_list:
        bas_mat[[row_index,i]]=bas_mat[[i,row_index]]
        row_index+=1
    for i in free_list:
        bas_mat[i,i]=1
    return A,b,bas_mat,free_list

if __name__ == "__main__":
    #Arg generation
    K=PolynomialRing(F,n*4,'x')
    x=K.gens()
    #Secret key
    key=vector([F(1) for _ in range(n)])
    key[0]=0
    key[1]=0
    key[2]=0
    #Roots of equations
    root=[]
    root.extend(key)
    #Linear polynomials
    linear_poly=[]
    #Numbar of variables
    var_cnt=0
    #Variables
    key_v=vector(x[0:n])
    var_cnt+=n
    #Input differences
    deta=vector([F(0) for _ in range(n)])
    deta[n-2]=F(1)
    deta[n-1]=F(1)
    
    #Arg for LowMC
    encrypt_mat_array,keygen_mat_array,round_constant_array,plain0=load_affine_layers(F,"./Affine_layers.txt")
    #The first plaintext and ciphertext pair
#     plain0=vector([F(0) for _ in range(n)])
#     for i in range(n):
#         plain0[i]=F.random_element()
    cipher0=LowMC(plain0,key,encrypt_mat_array,keygen_mat_array,round_constant_array,0,r)
    print("The first plain_cipher pair obtained")
    #The second plaintext and ciphertext pair
    plain1=plain0+deta
    cipher1=LowMC(plain1,key,encrypt_mat_array,keygen_mat_array,round_constant_array,0,r)
    print("The second plain_cipher pair obtained")
    Sbox_in_diff,Sbox_out_diff,active_cnt_list=LowMC_diff_path(plain0,plain1,key,encrypt_mat_array,keygen_mat_array,round_constant_array,0,r)
                
    #Internal states
    Sbox_in,Sbox_out=LowMC_state_path(plain0,key,encrypt_mat_array,keygen_mat_array,round_constant_array,0,r)
    #Generate equations
    state_out=cipher0
    state_out+=keygen_mat_array[r-1]*key_v
    state_out=affine_layer_inverse(encrypt_mat_array[r-1],round_constant_array[r-1],state_out)
    for i in range(r-1,0,-1):
        print("Equations of round_{} obtained".format(i))
        #Input state
        state_in=[]
        active_cnt=0
        for j in range(n//3):
            #Output difference
            Sbox_out_diff_index=vector_to_index(vector(Sbox_out_diff[i][j*3:j*3+3]))
            Sbox_in_diff_index=vector_to_index(vector(Sbox_in_diff[i][j*3:j*3+3]))
            #Inactive sboxs,represent input state with three variables
            if Sbox_out_diff_index==0 and Sbox_in_diff_index==0:
                if i==1:
                    continue
                state_in.extend(x[(var_cnt):(var_cnt+3)])
                #Add values of internal variables
                root.extend(Sbox_in[i][j*3:j*3+3])
                var_cnt+=3
            else:#Active sboxs,add two linear equations
                active_cnt+=1
                coe0=equ_t[Sbox_in_diff_index][Sbox_out_diff_index][0][0]
                coe1=equ_t[Sbox_in_diff_index][Sbox_out_diff_index][0][1]
                state=vector(state_out[j*3:j*3+3])
                constant=equ_t[Sbox_in_diff_index][Sbox_out_diff_index][1]
                linear_poly.extend([coe0*state+constant[0],coe1*state+constant[1]])
                #Calculate input state
                coe0=state_in_t[Sbox_in_diff_index][Sbox_out_diff_index][0][0]
                coe1=state_in_t[Sbox_in_diff_index][Sbox_out_diff_index][0][1]
                coe2=state_in_t[Sbox_in_diff_index][Sbox_out_diff_index][0][2]
                constant=state_in_t[Sbox_in_diff_index][Sbox_out_diff_index][1]
                state_in.extend([coe0*state+constant[0],coe1*state+constant[1],coe2*state+constant[2]])
        print("Round_{}:number of active S-boxes={}".format(i,active_cnt))
        #Calculate output state
        if i==1:
            continue
        state_out=vector(state_in)+keygen_mat_array[i-1]*key_v
        state_out=affine_layer_inverse(encrypt_mat_array[i-1],round_constant_array[i-1],state_out)
    root.extend([0 for _ in range(n*4-var_cnt)])
    test_poly(linear_poly,root)
    
    #Solve equations
    #Get coefficient matrix of equations
    mat,vec=poly_to_affine_mat(F,x,var_cnt,linear_poly)
    #Gaussian_elimination
    mat,vec,bas_mat,free_var_list=gaussian_elimination(F,mat,vec,mat.ncols())
    var_vec=[]
    var_vec.extend(x[0:var_cnt])
    var_vec.append(1)
    var_vec=vector(var_vec)
    result=bas_mat*var_vec
    test_result(result,root)
    for i in range(len(result)):
        print("x{}={}".format(i,result[i]))