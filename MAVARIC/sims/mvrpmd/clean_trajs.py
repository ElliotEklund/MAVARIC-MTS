import numpy as np

def clean_me(tag,num_beads):
  trajs_file = "Output/Trajectories/" + tag
  broken_file = "Output/" + "broken"

  v_traj = np.loadtxt(trajs_file)
  v_broke = np.loadtxt(broken_file)

  num_trajs = len(v_traj)/num_beads
  num_broke = len(v_broke)
  v_clean = np.zeros((num_trajs-num_broke)*num_beads)

  count_broke = 0
  count_clean = 0
  traj = 0

  while count_broke < num_broke:
    if(traj == v_broke[count_broke]):
      count_broke += 1
    else:
      for bead in range(num_beads):
        v_clean[count_clean] = v_traj[traj*num_beads + bead]
        count_clean += 1
    traj +=1

  while traj < num_trajs:
    for bead in range(num_beads):
      v_clean[count_clean] = v_traj[traj*num_beads + bead]
      count_clean += 1
    traj +=1

  return v_clean


def clean_me_elec(tag,num_beads,num_states):
  trajs_file = "Output/Trajectories/" + tag
  broken_file = "Output/" + "broken"

  v_traj = np.loadtxt(trajs_file)
  v_broke = np.loadtxt(broken_file)

  num_trajs = len(v_traj)/(num_beads*num_states)
  num_broke = len(v_broke)
  v_clean = np.zeros((num_trajs-num_broke)*num_beads*num_states)


  count_broke = 0
  count_clean = 0
  traj = 0
  while count_broke < num_broke:
    if(traj == v_broke[count_broke]):
      count_broke += 1
    else:
      for bead in range(num_beads):
        for state in range(num_states):
          v_clean[count_clean] = v_traj[traj*num_beads*num_states + bead*num_states + state]
          count_clean += 1
    traj +=1


  while traj < num_trajs:
    for bead in range(num_beads):
      for state in range(num_states):
        v_clean[count_clean] = v_traj[traj*num_beads*num_states + bead*num_states + state]
        count_clean += 1
    traj +=1


  return v_clean

#main
nuc_beads = 6
elec_beads = 1
num_states = 2

tag_Q = "Q"
tag_P = "P"
tag_x = "xElec"
tag_p = "pElec"

Q_clean = clean_me(tag_Q,nuc_beads)
P_clean = clean_me(tag_P,nuc_beads)
x_clean = clean_me_elec(tag_x,elec_beads,num_states)
p_clean = clean_me_elec(tag_p,elec_beads,num_states)

np.savetxt("Output/Trajectories/" + tag_Q + "clean",Q_clean)
np.savetxt("Output/Trajectories/" + tag_P + "clean",P_clean)
np.savetxt("Output/Trajectories/" + tag_x + "clean",x_clean)
np.savetxt("Output/Trajectories/" + tag_p + "clean",p_clean)

num_trajs_clean = len(Q_clean)/nuc_beads

print("Number of clean trajectories:" + num_trajs_clean)
