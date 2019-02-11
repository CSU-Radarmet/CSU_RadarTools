## CSU_RadarTools Release Instructions

To release a new version of `csu_radartools` on PyPI:

1. `git checkout master`

2. `git fetch upstream && get merge upstream/master`

3. `git clean -xfd`

4. Ensure CHANGELOG is up to date

5. Update `_version.py` (set release version, remove 'dev0')

6. `git add . && git commit -m 'Release X.X.X'`

7. `conda update pip wheel setuptools twine numpy cython` or `pip install -U pip` then `pip install -U wheel setuptools twine numpy cython`

8. `python setup.py sdist bdist_wheel`

9. `twine upload dist/*`

10. `git tag -a vX.X.X -m 'Release X.X.X'`

11. Update `_version.py` (add 'dev0' and increment minor)

12. `git add . && git commit -m 'Increment to dev version'`

13. `git checkout master`

14. `git push upstream master`

15. `git push upstream --tags`

16. Publish release annoucement to Github Releases page
