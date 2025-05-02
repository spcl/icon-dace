------------------------------------------------------
DWD Setup for BuildBot (on SX-AURORA)
------------------------------------------------------

[2017-07-31] : F. Prill, DWD
[2020-05-22] : D. Reinert, DWD: Adapted to SX-AURORA
[2020-09-21] : M Jacob, DWD: Ported most experiments to make_runscripts + make_target_runscript mechanism

------------------------------------------------------
DWD Setup for BuildBot
    Running the DWD Tests
    Running the DWD Tests manually
    Creating New DWD Tests
        Nomenclature
        Queue Restrictions for User "buildbot"
        Registering Tests in the "dwd" Experiment List
        Test Data and Output Directory
        Signalling "SUCCESS" or "FAILURE"
    Admin Information
        User "buildbot"
        Access Rights
        Directory Structure for BuildBot at DWD
------------------------------------------------------




h1. DWD Setup for BuildBot

The "dwd" experiment list contains a set of medium-size test cases which

* are executed on the DWD machine only (all other builders inactive),
* Most run scripts can be generated with make_target_runscript (Either through BuildBot or the make_runscripts front end)
* A few experiments are still plain PBSpro run scripts and do not use the usual BuildBot macros.


h2. Running the DWD Tests

Select the following parameters on the Buildbot web site (https://buildbot.dkrz.de/builders):

| Repository: | git@gitlab.dkrz.de:icon-nwp.git |
| Properties: | Name: "list" ; value: "dwd"     |


h2. Running the DWD Tests manually
*.run scripts in this directory can be submitted directly witch qsub.
For other tests see ../checksuite.nwp/ and the README.txt therein.


h2. Creating New DWD Tests

The DWD experiment list "dwd" uses the "plain run script mechanism" that has been implemented by Leonidas and Ralf.
When importing existing scripts, only a few adjustments are necessary.


h3. Nomenclature

Plain run scripts are located in <pre>run/checksuite.rcnl.dwd.de</pre>.

# It is necessary to rename the test scripts according to the following scheme: <pre>exp.<test name>.run</pre>
# It is necessary to set the stdout log output in the PBS header according to the following scheme: 
<pre>
#PBS -j oe
#PBS -o LOG.exp.<test name>.run.run.o 
</pre>

The scripts are started from the subdirectory @run@ with the help of a soft link to the @run/checksuite.rcnl.dwd.de@ subdirectory.
Therefore the environment variable @${PBS_O_WORKDIR}@ points to the @run@ subdirectory. 
This must be taken into account when using @${PBS_O_WORKDIR}@ together with relative paths.

Other run scripts templates are located in <pre>run/checksuite.nwp</pre>.
Those templates are converted into full run scripts using <pre>make_target_runscript</pre>.
Also see <pre>../checksuite.nwp/README.txt</pre>


h3. Queue Restrictions for User "buildbot"

* so far there are no queue restrictions on @sx_norm@ for the user buildbot

The current settings can be viewed by typing
<pre>
qstat -Q -f <queue name>
</pre>


h3. Registering Tests in the "dwd" Experiment List

New tests can be registered by extending the script
scripts/buildbot_scripts/create_list_dwd
as follows.
For template job scripts:
<pre>
$addexp checksuite.nwp/<template name> --builders DWD_nec --list $listname
</pre>
Note that <template name> is the file name, including the prefix "nwpexp." but without the suffix ".run".

For plain job scripts:
<pre>
$addexp checksuite.rcnl.dwd.de/<test name> --builders DWD_nec --list $listname
</pre>
Note that <test name> is the full name, including the prefix "exp." and the suffix ".run".



h3. Test Data and Output Directory

Test data (initial state, boundary data etc.) may be loaded from the directory <pre>/hpc/rwork0/routfor/test/icon/Checksuite_data</pre>.
Please ask T. Hanisch, e.g., for setting up a new subdirectory for your test case.

Output directories for the experiments may be created in @/hpc/gtmp/buildbot@ (available via the environment variable @$TMPDIR@) where they are automatically deleted.
The BuildBot system does not take care of cleaning your experiment directory before or after each run. You may take precautions in your own run scripts for this.


h3. Signalling "SUCCESS" or "FAILURE"

Run scripts signal "SUCCESS" or "FAILURE" ("green" or "red") with the help of the script error status:
Ending the script with
<pre>
 exit 0
</pre>
(or omitting the explicit @exit@ call) means that the script has run successfully.
All other @exit@ values mean "FAILURE".


h2. Admin Information


h3. User "buildbot"

A user "buildbot" is required. This user needs shell login, otherwise it is not possible to execute run scripts (though it is possible to submit them to the batch queue).
At DWD the user "buildbot" has the unusual home directory @/hpc/sw/buildbot/home@. 
Make sure that the run scripts require a proper setting of the environment variable @$TMPDIR@.


h3. Access Rights

For security reasons, the user "buildbot" is not part of the group FE1 and therefore does not have read or write access to these parts of the file system.
Only the read access to the necessary input data is granted by setting access control lists (ACLs) for the subdirectory
<pre>
/hpc/rwork0/routfor/test/icon/Checksuite_data
</pre>
by typing
<pre>
setfacl -d -R -m g:buildbot:rX /lustre2/rwork0/routfor/test/icon/Checksuite_data
setfacl  -R -m g:buildbot:rX /lustre2/rwork0/routfor/test/icon/Checksuite_data

getfacl /lustre2/rwork0/routfor/test/icon/Checksuite_data
</pre>


h3. Directory Structure for BuildBot at DWD

The BuildBot runs reside in the directory
<pre>
/hpc/uwork/buildbot/icon/nec/default
</pre>
