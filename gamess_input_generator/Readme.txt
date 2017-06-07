I have included some notes on common codes used in the Zimmerman Group, feel free to add additional
information and update index as needed!

Index                   - Line 4
Working Directory Files - Line 12
Basis Set Configuration - Line 25
Header Configuration    - Line 91
  Formatting            - Line 113
  Keywords              - Line 120
Operating the Program   - Line 128
Relevant Gamess options - Line
Other Thoughts          - 

############## Working Directory Files ################

First a foremost, make sure that in the working directory you have the
following files:

gamesinp.qsh
header
make_games_input.py
runall_py

As well as the following folder:
basis_sets

############### Basis Set Configuration ###################

In the basis_sets folder you should have your basis sets formatted as follows
(see below for some help on this, do not include lines with hashtags
preceeding them):

#For basis sets (ex. C 6-311++G**):#

S   6
  1   4563.2400000              0.00196665       
  2    682.0240000              0.0152306        
  3    154.9730000              0.0761269        
  4     44.4553000              0.2608010        
  5     13.0290000              0.6164620        
  6      1.8277300              0.2210060        
L   3
  1     20.9642000              0.1146600              0.0402487        
  2      4.8033100              0.9199990              0.2375940        
  3      1.4593300             -0.00303068             0.8158540        
L   1
  1      0.4834560              1.0000000              1.0000000        
L   1
  1      0.1455850              1.0000000              1.0000000        
L   1
  1      0.0438000              1.0000000              1.0000000        
D   1
  1      0.6260000              1.0000000 

#For ECP:#

  PD-ECP GEN     28     3   
  5      ----- f-ul potential     -----
     -0.0563177        0    598.3336444    
    -20.1288036        1    162.4298290    
   -105.8197923        2     51.5714771    
    -42.5733345        2     16.4888260    
     -3.6165086        2      5.8287656    
  5      ----- s-ul potential     -----
      3.0003651        0     73.3806304    
     32.4350093        1     14.7550438    
    459.0830383        2     17.8350204            
   -868.0629029        2     12.7111477            
    514.4726098        2      9.3292063            
  5      ----- p-ul potential     -----
      4.9593099        0     55.6689877            
     21.1711029        1     64.2337771            
    605.0560092        2     17.6254952            
   -726.9641846        2     11.9058155            
    396.3274883        2      8.5100832            
  5      ----- d-ul potential     -----
      3.0508745        0     49.9994728            
     22.2506580        1     39.7477547            
    674.8357698        2     11.4321366            
  -1040.8554048        2      9.1790080            
    505.9375147        2      7.5624429            

Basis sets can be obtained from https://bse.pnl.gov/bse/portal (basis set
exchange).  This website has a pull down menu called "Format", under this menu
select the GAMESS-US output format and then proceed to select relevant
atoms/basis set.  
*****Minor modifications to this GAMESS-US default format is
neccessary to match the above format, the program WILL NOT WORK unless you
match the format shown******

######## HEADER CONFIGURATION ###############

Detailed documentation on different keywords for the Gamess header can be
found at http://www.msg.ameslab.gov/gamess/GAMESS_Manual/input.pdf, or in the
Zimmerman group dropbox folder.  If you do not have access to the group
dropbox folder please see another group member for help.

An example header file is included in the SVN repository and looks as follows:

##
 $CONTRL DFTTYP=B3LYP RUNTYP=ENERGY NOSYM=1
   UNITS=ANGS PP=READ MAXIT=200 
   ICHARG=0 $END
 $SYSTEM MWORDS=500 $END
 $SCF DIIS=.TRUE. $END
! $PCM SOLVNT=ACETNTRL SMD=.TRUE. $END
 $PCM SOLVNT=DCM SMD=.TRUE. $END
 $DATA

C1
##

A couple notes on FORMATING:
The spaces preceding the individual lines of the header file are 100%
neccessary for functionality of Gamess.  I don't know why, something to do
with fortran programming language. Adding a ! prior to a line comments it out.
Make sure you do not add extra blank lines to the start or end of the header
file. 

A couple notes on KEYWORDS:
DFTTYP=B3LYP  ## This sets the level of theory
RUNTYP=ENERGY ## sets job type
PP=READ       ## indicates that the job is ECP !!!!! PLEASE MANUALLY REMOVE
                 THIS IF YOU ARE NOT RUNNING ECP !!!!
ICHARG=0      ## sets the charge of the overall system
SOLVNT=???    ## sets solvent of your choice see keywords section for more

########## Operating the Program ######### 

The program is fairly straight forward:
1) Import XYZ's into working directory with all needed files for this program
2) Double check that the basis_set folder contains proper basis set
3) Make sure header file contains keywords for YOUR needs
4) If you are running ECP runall_py is fine
	IF NOT change the line:
   python make_games_input.py -ECP 1 $item 
        TO
   python make_games_input.py -ECP 0 $item 
5) ./runall_py  # this will generate input files from your XYZ files using the
                  basis set and header that you have placed in the working 
                  directory
6) For more details on operating the program, please open the file
make_gamess_input.py in a text editor (vim?) and read the first few lines of
this file.

########## Relevant Gamess Options ######## 
Please add newly discovered or useful Gamess options here for all other
options please refer to gamess documentation for more specifics!

 


http://myweb.liu.edu/~nmatsuna/gamess/input/PCM.html
