'''Script that updates the OpenMM imports in Duck scripts once it has been installed in the conda environment
args:
  - rootDir: environment directory containing the scripts (it is passed from the Scipion3 installer)

'''

import sys, os

endDir = 'duck/steps'
OLD_STR, NEWSTR = 'Platform_getPlatformByName', 'Platform.getPlatformByName'
REPLACE_DIC = {f'{endDir}/equlibrate.py': [OLD_STR, NEWSTR],
               f'{endDir}/normal_md.py': [OLD_STR, NEWSTR],
               f'{endDir}/steered_md.py': [OLD_STR, NEWSTR], }


def findFullPath(rootDir, targetSuffix):
  for dirpath, dirnames, filenames in os.walk(rootDir):
    if 'openduck' in dirpath:
      for filename in filenames:
        fullFilePath = os.path.join(dirpath, filename)
        if fullFilePath.endswith(targetSuffix):
          return os.path.abspath(fullFilePath)
  return None

def replaceInFile(file, inStr, outStr):
  with open(file) as f:
    text = f.read()

  nText = text.replace(inStr, outStr)
  with open(file, 'w') as f:
    f.write(nText)

if __name__ == "__main__":
  rootDir = sys.argv[1]
  repComds = []
  for suffixPath, reps in REPLACE_DIC.items():
    scriptPath = findFullPath(rootDir, suffixPath)
    if scriptPath:
        replaceInFile(scriptPath, reps[0], reps[1])
    