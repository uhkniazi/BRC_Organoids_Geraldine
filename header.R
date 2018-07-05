# File: header.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: global variables
# Date: 05/07/2018


## variables
g_pid = 13
g_did = 28
gcswd = getwd()
gcRemoteDir = "/run/user/1000/gvfs/sftp:host=10.202.64.29,user=k1625253/users/k1625253/brc_scratch/Data/ProjectsData/BRC_Organoids_Geraldine/"

p.old = par()

###### utility functions

f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}