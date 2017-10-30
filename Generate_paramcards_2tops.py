
# Generate param cards for EFT runs for M_2tops for Lambda = 1,2,3,4

# RUN1
run1_c1 = -1.75
run1_c2 = -1.76
#rest 0

# RUN2
run2_c1 = 1.62
run2_c2 = 1.61
#rest 0

# RUN3
run3_c3 = 4.90
run3_c4 = 2.94
#rest 0

# RUN4
run4_c3 = 4.90
run4_c5 = 5.46
#rest 0

# RUN5
run5_c1 = 1.62
run5_c3 = 4.90
#rest 0



for LAMBDA in range(1,5):
    LAMBDA2 = LAMBDA*LAMBDA

    filename = [""]*5
    for n_run in range(1,6):
        filename[n_run-1] = "param_card_LAMBDA"+str(LAMBDA)+"_run"+str(n_run)+".dat"

    ##### RUN1 ######
    fin = open("param_card.dat")
    fout_run1 = open(filename[0],"wt") 

    # Copy first part unchanged of param card
    for i in range(20): # WORKS!
        fout_run1.write(fin.readline())

    # write correct run numbers for c_i
    fout_run1.write("    1 "+str(run1_c1*LAMBDA2)+" # CRL2 \n")
    fout_run1.write("    2 "+str(run1_c2*LAMBDA2)+" # C1LL2 \n")
    fout_run1.write("    3 "+str(0)+" # C8LL2 \n")
    fout_run1.write("    4 "+str(0)+" # C1BL2 \n")
    fout_run1.write("    5 "+str(0)+" # C8BL2 \n")

    # Skip over 5 lines
    for i in range(5):
        fin.readline()

    # Copy rest of file
    for i in range(100): # WORKS!
        fout_run1.write(fin.readline())

    # Close both files
    fout_run1.close()
    fin.close()


    ##### RUN2 ######
    fin = open("param_card.dat")
    fout_run2 = open(filename[1],"wt") 

    # Copy first part unchanged of param card
    for i in range(20): # WORKS!
        fout_run2.write(fin.readline())

    # write correct run numbers for c_i
    fout_run2.write("    1 "+str(run2_c1*LAMBDA2)+" # CRL2 \n")
    fout_run2.write("    2 "+str(run2_c2*LAMBDA2)+" # C1LL2 \n")
    fout_run2.write("    3 "+str(0)+" # C8LL2 \n")
    fout_run2.write("    4 "+str(0)+" # C1BL2 \n")
    fout_run2.write("    5 "+str(0)+" # C8BL2 \n")

    # Skip over 5 lines
    for i in range(5):
        fin.readline()

    for i in range(100): # WORKS!
        fout_run2.write(fin.readline())
        
    fout_run2.close()
    fin.close()



        ##### RUN 3 ######
    fin = open("param_card.dat")
    fout_run3 = open(filename[2],"wt") 

    # Copy first part unchanged of param card
    for i in range(20): # WORKS!
        fout_run3.write(fin.readline())

    # write correct run numbers for c_i
    fout_run3.write("    1 "+str(0)+" # CRL2 \n")
    fout_run3.write("    2 "+str(0)+" # C1LL2 \n")
    fout_run3.write("    3 "+str(run3_c3*LAMBDA2)+" # C8LL2 \n")
    fout_run3.write("    4 "+str(run3_c4*LAMBDA2)+" # C1BL2 \n")
    fout_run3.write("    5 "+str(0)+" # C8BL2 \n")

    # Skip over 5 lines
    for i in range(5):
        fin.readline()

    for i in range(100): # WORKS!
        fout_run3.write(fin.readline())
        
    fout_run3.close()
    fin.close()




        ##### RUN4 ######
    fin = open("param_card.dat")
    fout_run4 = open(filename[3],"wt") 

    # Copy first part unchanged of param card
    for i in range(20): # WORKS!
        fout_run4.write(fin.readline())

    # write correct run numbers for c_i
    fout_run4.write("    1 "+str(0)+" # CRL2 \n")
    fout_run4.write("    2 "+str(0)+" # C1LL2 \n")
    fout_run4.write("    3 "+str(run4_c3*LAMBDA2)+" # C8LL2 \n")
    fout_run4.write("    4 "+str(0)+" # C1BL2 \n")
    fout_run4.write("    5 "+str(run4_c5*LAMBDA2)+" # C8BL2 \n")

    # Skip over 5 lines
    for i in range(5):
        fin.readline()

    for i in range(100): # WORKS!
        fout_run4.write(fin.readline())
        
    fout_run4.close()
    fin.close()




        ##### RUN5 ######
    fin = open("param_card.dat")
    fout_run5 = open(filename[4],"wt") 

    # Copy first part unchanged of param card
    for i in range(20): # WORKS!
        fout_run5.write(fin.readline())

    # write correct run numbers for c_i
    fout_run5.write("    1 "+str(run5_c1*LAMBDA2)+" # CRL2 \n")
    fout_run5.write("    2 "+str(0)+" # C1LL2 \n")
    fout_run5.write("    3 "+str(run5_c3*LAMBDA2)+" # C8LL2 \n")
    fout_run5.write("    4 "+str(0)+" # C1BL2 \n")
    fout_run5.write("    5 "+str(0)+" # C8BL2 \n")

    # Skip over 5 lines
    for i in range(5):
        fin.readline()

    for i in range(100): # WORKS!
        fout_run5.write(fin.readline())
        
    fout_run5.close()
    fin.close()



    


        
