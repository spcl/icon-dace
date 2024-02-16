import subprocess
import unittest

from engine.run_ensemble import append_job, finalize_jobs


class TestRunEnsemble(unittest.TestCase):
    def testJobSuccessfulSerial(self):
        job = "test 1 = 1".split()
        append_job(job=job, job_list=[], parallel=False)

    def testJobSuccessfulParallel(self):
        job = "test 1 = 1".split()
        job_list = []
        append_job(job=job, job_list=job_list, parallel=True)
        finalize_jobs(job_list=job_list, dry=False, parallel=True)

    def testJobFailedSerial(self):
        job = "test 1 = 2".split()
        with self.assertRaises(subprocess.CalledProcessError):
            append_job(job=job, job_list=[], parallel=False)

    def testJobFailedParallel(self):
        job = "test 1 = 2".split()
        job_list = []
        append_job(job=job, job_list=job_list, parallel=True)
        with self.assertRaises(subprocess.CalledProcessError):
            finalize_jobs(job_list=job_list, dry=False, parallel=True)


if __name__ == "__main__":
    unittest.main()
