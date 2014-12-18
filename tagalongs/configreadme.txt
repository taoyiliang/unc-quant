If you are using an OpenSSH compatible client, you can place the following
file in your ~/.ssh/ directory.
<config>

To utilize this file you will first need change the userid in the config file
from USERNAME to your userid (Line numbers:8,9, and 15). Next to login to
hpclogin.inel.gov you only type hpclogin since the script will add the domain.
Once a connection has been established to hpclogin leave that session active.
You may then start additional SSH sessions directly to the machines within the
HPC enclave and they will be responsive and will not require you to enter your
RSA token additional times.
