# Author: Tanya Mehta, Parkland High School
# Cell Membranes represnted as 2 Dimensional network with 100 beads (10x10) sheets. Cell Membrane is moving, Virus Membrane is fixed
#Starting Position for Cell membrane at z=0
import pandas as pd
import math
import numpy as np
import Graph_def
from xlwt import Workbook 
#Open Workbook
wb = Workbook() 
#*************************************************************************************
# add_sheet is used to create sheet in the excel output file 
sheet1 = wb.add_sheet('CSI_Time')
sheet2 = wb.add_sheet('Output-CSI_AttachTime')
#Read Excel input File
#Original Position of the 100 cell nodes.
node_position_original = pd.read_excel("membrane_100.xlsx", sheet_name='Node_Position_Original')
# Connectivity of the nodes via springs and spring constants *x-direction, y-direction and diagonals)
spring = pd.read_excel("membrane_100.xlsx", sheet_name='Springs')
#Original Position of the 100 cell nodes after 1 % perturbation in x and y directions.
node_position_perturbed = pd.read_excel("membrane_100.xlsx", sheet_name='Node_Position_Perturbed')
# Membrane CG - Node and Spring Configuration - 10 rows x 10 columns. Number of nodes = 100, Number of springs = 40
num_nodes = len(node_position_original)
num_springs = len(spring)
columns = node_position_original['x'].max()+1
rows = node_position_original['y'].max()+1 
#length between the nodes 
Len_Node = 4.445e-9
TA = [[[0 for z in range(2)] for x in range(10)]for y in range(10)]
write_index=0
#Loop to study Temperature(K): (280,290,300,310,320) and Membrane Tension Force (pN): (0,50, 100, 150, 200, 250)
for o in range(0,5):
    for q in range(0,6):
# Simulation Parameters********************************** 
        timestep =  5e-13
        #Number of time steps
        num_iterations = 40000
        time_counter = 0.0
        #Resistance of the medium
        eta = 2.5e-13
        Epsilon=16e-21
        #Virus definitions
        #Shape is sphere
        #Radius of the virus
        R_V = 1e-7
        #Radius of the cell and virus beads
        R_CBead = 2.5e-9
        R_VBead = 1e-9
        #distance parameter in the potential function.
        Sigma=2.5e-9
        #Initial separation between virus and cell membrane
        Virus_sep = 6e-9
        #minimum separation distance
        ProximityMax = num_nodes*np.power(2,(1/6))*Sigma
        #dimension of the virus membrane - 1 bead representation of the virus
        Dim_V = 1
# Constants********************************** 
        kB = 1.38e-23 # Boltzman constant
        #Induced Tension
        Tension = (0+q*5e-11)
        #Temperature
        Temp = (280+10*o)
        #Boundary Condition- Induced lateral tension in x and y directions on these nodes
        fix_nodes_x = [0,9,10,19,20,29,30,39,40,49,50,59,60,69,70,79,80,89,90,99]
        fix_nodes_y = [0,90,1,91,2,92,3,93,4,94,5,95,6,96,7,97,8,98,9,99]
        #Initialization
        # s_0 is the original prositions of 100 nodes (100,2)
        s_0 = [[0 for x in range(3)] for y in range(num_nodes)]
        #Assign S_O the original node positions
        for i in range(0,num_nodes):
            s_0[i][0] = node_position_original.x[i]*Len_Node
            s_0[i][1] = node_position_original.y[i]*Len_Node
            s_0[i][2] = node_position_original.z[i]*Len_Node
        # s_p is the perturbed prositions of 25 nodes (100,2)
        s_p = [[0 for x in range(3)] for y in range(num_nodes)]
        #Position for 100 nodes at num_iterations time steps
        Position = [[[0 for x in range(3)] for y in range(num_nodes)] for z in range(num_iterations)]
        #Read or Calculate the Perturbed Node Positions and assign the original Position
        for i in range(0,num_nodes):
            s_p[i][0] = node_position_perturbed.x[i]*Len_Node
            s_p[i][1] = node_position_perturbed.y[i]*Len_Node
            s_p[i][2] = node_position_perturbed.z[i]*Len_Node
            Position[0][i][0] = s_p[i][0]
            Position[0][i][1] = s_p[i][1]
            Position[0][i][2] = s_p[i][2]
        #For Every spring the connecting nodes
        Nodes = [[0 for x in range(2)] for y in range(num_springs)]
        #Initialization of Arrays
        # Temp array for calculatining distance between coordinates
        dist = [[0 for x in range(3)] for y in range(num_springs)]
        # Original unperturbed length of springs Each os the spring is len_node of length
        Len_Orig = [0 for x in range(num_springs)]
        for i in range(0,num_springs):
            Len_Orig[i] = spring.Length[i]*Len_Node
        # Length of the spring at each time step.
        Len_Pert = [0 for x in range(num_springs)]
        theta = [0 for x in range(num_springs)]
        phi = [0 for x in range(num_springs)]
        # x_pos, y_pos, z_pos are temporary arrays for node coordinates to create different graphs 
        x_pos = [0.0 for x in range(num_nodes)]
        y_pos = [0.0 for x in range(num_nodes)]
        z_pos = [0.0 for x in range(num_nodes)]
        # x_direct, y_direct, z_direct are temporary arrays for spring force vectors to create different graphs 
        x_direct = [0 for x in range(num_nodes)]
        y_direct = [0 for x in range(num_nodes)]
        z_direct = [0 for x in range(num_nodes)]
        #color and nodenumber are temporary variables.
        color = [0 for x in range(num_nodes)]
        Proximity = [0 for x in range(num_iterations)]
        nodenumber = [0 for x in range(num_nodes)]
        #Spring Force for each of the springs
        F_U = [[0 for x in range(3)] for y in range(num_springs)]
        #unbalanced force per node
        F_Node = [[0 for x in range(3)] for y in range(num_nodes)]
        #Random Force on each of the nodes
        F_Random = [[0 for x in range(3)] for y in range(num_nodes)]
        #Virus Interaction Force on Cell
        F_Cell = [[0 for x in range(3)] for y in range(num_nodes)]
        #Unbalanced Force for each nodes at num_iterations time steps
        Unbalanced_Force = [[[0 for x in range(3)] for y in range(num_nodes)] for z in range(num_iterations)]
        # Calculate Force due to temperature
        R = np.sqrt(6*eta*kB*Temp/timestep)   
        #Number of beads on the virus
        num_nodes_V = Dim_V*Dim_V
        #Virus Interaction Force on Virus due to cell beads
        F_V = [[0 for x in range(3)] for y in range(num_nodes_V)]
        # Distance vector between two beads. 
        dist1 = [[0 for x in range(3)] for y in range(num_nodes_V)]
        dist2 = [0 for x in range(num_nodes)]
        #Theta and phi angles for unit vetor between two beads
        thetaV = [0 for x in range(num_nodes_V)]
        phiV = [0 for x in range(num_nodes_V)]
        #magnitude of cell-virus interaction force between each bead and Virus Membrane.
        FVC_List = [[0 for x in range(num_nodes)] for y in range(num_nodes_V)]
        #Initial placent of virus relative to cell. 
        dist_V = [0.0 for x in range(3)]
        dist_V[0] = ((np.power(num_nodes,0.5)-1)/2)*Len_Node
        dist_V[1] = ((np.power(num_nodes,0.5)-1)/2)*Len_Node
        #Distance of Virus from Membrane (Perpendicular distance from center to the membrane)
        dist_V[2] = Virus_sep
        #Location of Virus beads at different time steps(It is not changing in this simulation)
        Location_Node_V = [[[0 for x in range(3)] for y in range(num_nodes_V)] for z in range(num_iterations)]
        Separation = [[0 for x in range(num_nodes)] for y in range(num_iterations)]
        Separation_Min = [0 for x in range(num_iterations)]
        Counter = [0 for x in range(num_iterations)]
        Sep = [0 for x in range(num_iterations)]
        # x_pos_V, y_pos_V, z_pos_V are temporary arrays for virus node coordinates to create different graphs 
        x_pos_V = [0.0 for x in range(num_nodes_V)]
        y_pos_V = [0.0 for x in range(num_nodes_V)]
        z_pos_V = [0.0 for x in range(num_nodes_V)]
        # x_pos_VC, y_pos_VC, z_pos_VC are temporary arrays for (virus_cell) nodes to create different graphs
        x_pos_VC = [0.0 for x in range(num_nodes_V+num_nodes)]
        y_pos_VC = [0.0 for x in range(num_nodes_V+num_nodes)]
        z_pos_VC = [0.0 for x in range(num_nodes_V+num_nodes)]
        for i in range(0,num_nodes):
                for j in range(0,Dim_V):
                    Location_Node_V[0][0][0] = dist_V[0] #Virus positioned in the center of the membrane
                    Location_Node_V[0][0][1] = dist_V[1] # Virus positioned in the center of the membrane
                    Location_Node_V[0][0][2] = R_V+ dist_V[2] 
        for j in range(0,num_nodes):
            for i in range(0,num_nodes_V):
                k = i+j*num_nodes_V
                dist1[i][0] = (Position[0][j][0]-Location_Node_V[0][i][0])
                dist1[i][1] = (Position[0][j][1]-Location_Node_V[0][i][1])
                dist1[i][2] = (Position[0][j][2]-Location_Node_V[0][i][2])
                Separation[0][j] = (np.power(np.power(dist1[i][0],2.0)+np.power(dist1[i][1],2.0)+np.power(dist1[i][2],2.0),0.5)-R_V)
                Proximity[0]= Proximity[0] + (np.power(np.power(dist1[i][0],2.0)+np.power(dist1[i][1],2.0)+np.power(dist1[i][2],2.0),0.5)-R_V)
        # Start of time integration loop
        for l in range(1,num_iterations):
            time_counter = time_counter + timestep
        #Virus Interaction Force Vectors - Resetting at each iteration step
            F_Cell = [[0 for x in range(3)] for y in range(num_nodes)]
            F_V = [[0 for x in range(3)] for y in range(num_nodes_V)]
            FVC_List = [[0 for x in range(num_nodes)] for y in range(num_nodes_V)]
            for i in range(0,num_nodes):
        #Generation of random direction for the thermal fluctuation force vector
                R_theta = np.random.uniform(0,2*np.pi-1e-5)
        #Correcting for bias in Sperhical sampling
                v = np.random.uniform(0,0.999999999)
                R_phi = np.arccos(2*v-1)
        # Three Components of the Random Force
                F_Random[i][0] = R*np.cos(R_theta)*np.sin(R_phi)
                F_Random[i][1] = R*np.sin(R_theta)*np.sin(R_phi)
                F_Random[i][2] = R*np.cos(R_phi)
        #loop over each virus bead to calculate pairwise force
                for j in range(0,num_nodes_V):
                    dist1[j][0] = (s_p[i][0]-Location_Node_V[l-1][j][0])
                    dist1[j][1] = (s_p[i][1]-Location_Node_V[l-1][j][1])
                    dist1[j][2] = (s_p[i][2]-Location_Node_V[l-1][j][2])
                    distance = np.power(np.power(dist1[j][0],2.0)+np.power(dist1[j][1],2.0)+np.power(dist1[j][2],2.0),0.5)
        #Force calculation
                    FVC = 4*Epsilon*((12*np.power(Sigma,12)/np.power((distance-R_V),13))-(6*np.power(Sigma,6)/np.power((distance-R_V),7)))
                    FVC_List[j][i] = FVC
        #Three component vectors of force
                    F_Cell[i][0] = F_Cell[i][0] + FVC*(dist1[j][0]/distance)
                    F_Cell[i][1] = F_Cell[i][1] + FVC*(dist1[j][1]/distance)
                    F_Cell[i][2] = F_Cell[i][2] + FVC*(dist1[j][2]/distance)
            F_Node = [[0 for x in range(3)] for y in range(num_nodes)]
        #Calculate perturbations and force for each of the springs
            for i in range(0,num_springs):
                Nodes[i][0] = spring.Node1[i]  #assigning nodes for each spring 
                Nodes[i][1] = spring.Node2[i]
                j = Nodes[i][0] -1 # adjust for the index starting at zero for the first node
                k = Nodes[i][1] -1
        #Unperturbed Length of the spring
                dist[i][0] = (s_p[k][0]-s_p[j][0])
                dist[i][1] = (s_p[k][1]-s_p[j][1])
                dist[i][2] = (s_p[k][2]-s_p[j][2])
        #Perturbed Length of the spring
                Len_Pert[i] = np.power(np.power(dist[i][0],2.0)+np.power(dist[i][1],2.0)+np.power(dist[i][2],2.0),0.5)
        #Calculate three components of the spring Force on each of the springs  
                F_U[i][0] = -(dist[i][0]/Len_Pert[i])*spring.SpringConstant[i]*(Len_Pert[i]-Len_Orig[i])
                F_U[i][1] = -(dist[i][1]/Len_Pert[i])*spring.SpringConstant[i]*(Len_Pert[i]-Len_Orig[i])
                F_U[i][2] = -(dist[i][2]/Len_Pert[i])*spring.SpringConstant[i]*(Len_Pert[i]-Len_Orig[i])
        #Adding up the forces for all the springs on a given node
            for i in range(0,num_springs):
                j = Nodes[i][0] -1
                k = Nodes[i][1] -1
                F_Node[j][0] = F_Node[j][0]-F_U[i][0]
                F_Node[j][1] = F_Node[j][1]-F_U[i][1]
                F_Node[j][2] = F_Node[j][2]-F_U[i][2]
                F_Node[k][0] = F_Node[k][0]+F_U[i][0]
                F_Node[k][1] = F_Node[k][1]+F_U[i][1] 
                F_Node[k][2] = F_Node[k][2]+F_U[i][2]
        #Assigning Tension at the boundary nodes in x and y directions
            for i in range(0,20):
                F_Node[fix_nodes_x[i]][0] = F_Node[fix_nodes_x[i]][0] + np.power(-1,i+1)*Tension
                F_Node[fix_nodes_y[i]][1] = F_Node[fix_nodes_y[i]][1] + np.power(-1,i+1)*Tension
        #Storing the forces history 
            for j in range(0,num_nodes):
                Unbalanced_Force[l][j][0] = F_Node[j][0]
                Unbalanced_Force[l][j][1] = F_Node[j][1]
                Unbalanced_Force[l][j][2] = F_Node[j][2]
        # Update the position of the nodes (Euler Update)
            for j in range(0,num_nodes):
                Position[l][j][0]  = Position[l-1][j][0] + (F_Node[j][0]+F_Cell[j][0]+ F_Random[j][0])*timestep/eta 
                Position[l][j][1]  = Position[l-1][j][1] + (F_Node[j][1]+F_Cell[j][1]+ F_Random[j][1])*timestep/eta
                Position[l][j][2]  = Position[l-1][j][2] + (F_Node[j][2]+F_Cell[j][2]+ F_Random[j][2])*timestep/eta
        # Virus nodes ar not moving
            for j in range(0,num_nodes_V):
                Location_Node_V[l][j][0] = Location_Node_V[l-1][j][0]
                Location_Node_V[l][j][1] = Location_Node_V[l-1][j][1]
                Location_Node_V[l][j][2] = Location_Node_V[l-1][j][2]
        #Update the position of the perturbed nodes
            for j in range(0,num_nodes):
                s_p[j][0] = Position[l][j][0]
                s_p[j][1] = Position[l][j][1]
                s_p[j][2] = Position[l][j][2]
        #Graphical Output at every 500th iteration
            if ((l%500) == 0):
                for j in range(0,num_nodes):
                    nodenumber[j] = j+1
                    x_pos[j] = Position[l][j][0]
                    y_pos[j] = Position[l][j][1]
                    z_pos[j] = Position[l][j][2]
                Graph_def.UserPlot3DSurf(x_pos, y_pos, z_pos)
            for j in range(0,num_nodes):
                for i in range(0,num_nodes_V):
                    k = i+j*num_nodes_V
                    dist1[i][0] = (Position[l][j][0]-Location_Node_V[l][i][0])
                    dist1[i][1] = (Position[l][j][1]-Location_Node_V[l][i][1])
                    dist1[i][2] = (Position[l][j][2]-Location_Node_V[l][i][2])
                    Separation[l][j] = (np.power(np.power(dist1[i][0],2.0)+np.power(dist1[i][1],2.0)+np.power(dist1[i][2],2.0),0.5)-R_V)
                    Proximity[l]= Proximity[l] + (np.power(np.power(dist1[i][0],2.0)+np.power(dist1[i][1],2.0)+np.power(dist1[i][2],2.0),0.5)-R_V)
            Proximity[l]=  Proximity[l]/ProximityMax
            Separation_Min[l] = np.min(Separation[l])
        #End of time loop
        for j in range(0,num_nodes):
            x_pos[j] = Position[l][j][0]
            y_pos[j] = Position[l][j][1]
            z_pos[j] = Position[l][j][2]
            x_direct[j] = Unbalanced_Force[l][j][0]
            y_direct[j] = Unbalanced_Force[l][j][1]
            z_direct[j] = Unbalanced_Force[l][j][2]
            color[j] = math.sqrt(math.pow(x_direct[j],2)+math.pow(y_direct[j],2)+math.pow(z_direct[j],2))
        Graph_def.UserPlot3DSurf(x_pos, y_pos, z_pos)
        x_iter_tstep = [0.0 for x in range(num_iterations)]
        for i in range(0,num_iterations):
            x_iter_tstep[i]=timestep*(i)
        # Moving Average Window
        Window = 500
        Separation_Min_MA = [0.0 for x in range(num_iterations)]
        Graph_def.MA(Separation_Min_MA,Separation_Min,num_iterations,Window)
        ProximityMA = [0.0 for x in range(num_iterations)]
        Graph_def.MA(ProximityMA,Proximity,num_iterations,Window)
        SM_Average= np.average(Separation_Min_MA[(num_iterations-Window-10):(num_iterations-Window)])
        Proximity_Average= np.average(Proximity[(num_iterations-Window-10):(num_iterations-Window)])
        Attach_Index = 0
        Time_To_Attach = 0.0
        sheet1.write(write_index,0,Temp)
        sheet1.write(write_index,1,Tension)
        write_index = write_index +1
        for i in range(0,num_iterations):
            sheet1.write(write_index+i,0,x_iter_tstep[i])
            sheet1.write(write_index+i,1,ProximityMA[i])
            if ((ProximityMA[i] < 1.2) and (Attach_Index < 0.5)):
                Attach_Index = 1
                Time_To_Attach = x_iter_tstep[i]
        write_index = write_index + num_iterations
        wb.save('locationout_3D.xls')
        print(time_counter, Temp, Tension)
        Graph_def.UserPlot2Dline(x_iter_tstep,Sep,"Distance from Surface","Time (Seconds)","Distance from Surface (m)")
        Graph_def.UserPlot2Dline(x_iter_tstep,Proximity,"Cummulative Separation Index","Time (Seconds)","Proximity (m)")
        Graph_def.UserPlot2Dline(x_iter_tstep,Separation_Min,"Closest Point","Time (Seconds)","Closest Distance (m)")
        Graph_def.UserPlot2Dline(x_iter_tstep,Separation_Min_MA,"Moving Avg. of Closest Point","Time (Seconds)","Closest Distance (m)")
        Graph_def.UserPlot2Dline(x_iter_tstep,ProximityMA,"Moving Avg. of Cummulative Separation Index","Time (Seconds)","Proximity (m)")
        TA[o][q][0] = Time_To_Attach
        TA[o][q][1] = ProximityMA[num_iterations-Window-1]
#Excel Output to 'locationout_3D.xls'
for o in range(0,5):
    for q in range(0,6):
        k = (o)*6+q
        sheet2.write(k,0,280+o*10)
        sheet2.write(k,1,0+q*5e-11)
        sheet2.write(k,2,TA[o][q][0])
        sheet2.write(k,3,TA[o][q][1])
wb.save('locationout_3D.xls')
