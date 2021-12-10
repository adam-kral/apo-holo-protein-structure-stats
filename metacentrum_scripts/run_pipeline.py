"""

job.sh:

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/brno3-cerit/home/jenicek/test_directory # substitute username and path to to your real username and path

test -n "$SCRATCHDIR"

qsub -l select=1:ncpus=1:mem=4gb:scratch_local=100gb

echo $SCRATCHDIR
/scratch/melounova/job_9309132.meta-pbs.metacentrum.cz
"""
import string
from pathlib import Path

""" 
Quota/num jobs
(less than 10 minutes) jobs, we strongly recommend to run them in batches submitted as one job. To prevent PBS server 
glutting, there is a quota of 10 000 jobs (running or queuing) per user. 

Udelam treba 30 jobu? 5K-10K struktur na jednu, to sedi
# filter structures
# isoform budou asi dalsi joby.. (staci min cpu, to teda ale staci i na download, ale necham to tak)
# ani nemusim limitovat struktrutry (less than 10 minutes) jobs, we strongly recommend to run them in batches submitted as one job. To prevent PBS server glutting, there is a quota of 10 000 jobs (running or queuing) per user. 

Páry
# tam muze byt i vic jobů - muze jich být O(n^2), kde n je 
# to si vyzkousim lokalne, kolik jich bude


Storage
pouzij treba: (alzac)
brno12-cerit, unlimited quota

jeden z nich

ale zapisuje se vzdy do scratchdiru
a co pip python? Nainstalovat do jednoho? Asi ok, kdyz to je namapovany?
Spis ne, pokazdy je ten home dir taky jinej (asi podle jobu.. ale stejne se to posila po siti?) 

Asi ten pip install byl strasne pomaly, protoze se ty install soubory posilaly po siti!!!

/storage/brno11-elixir/home/baobab  # asi muj adresar na gridu, ale tech serveru je vic ne? Ma cenu ho pouzivat? Pro venv/
condu asi ne..

"""

"""
module add python/3.8.0-gcc-rab6t
python -m venv
activate..
python -m ensurepip --upgrade  # upgrade uplne nefunguval?
pip install --upgrade pip  # proto
pip install -r requirements.txt
# ted uz by melo vse fungovat
# spusťt skript s tim input souborem (cesta asi v argumentu jobu?)
qsub -v a=10, "var2='A,B'", c=20,  # var bude cesta k shardu (od datadir treba?)
# internet je strasne rychlej, asi bych mel mit min download threadů? Treba 10?

qsub -S python3 zkusit?

"""


""""""


class ShellTemplate(string.Template):
    delimiter = '><>'


def get_shell_template():
    with (Path(__file__).parent / 'generic_run_python_script.sh').open() as f:
        return ShellTemplate(f.read())
