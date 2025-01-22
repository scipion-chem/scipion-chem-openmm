'''Script that updates the OpenMM imports in Duck scripts once it has been installed in the conda environment
args:
  - rootDir: environment directory containing the scripts (it is passed from the Scipion3 installer)

'''

import sys, os

endDir = 'duck/steps'
REPLACE_DIC = {f'{endDir}/equlibrate.py': ['Platform_getPlatformByName', 'Platform.getPlatformByName'],
               f'{endDir}/normal_md.py': ['Platform_getPlatformByName', 'Platform.getPlatformByName'],
               f'{endDir}/steered_md.py': ['Platform_getPlatformByName', 'Platform.getPlatformByName'], }


def findFullPath(root_dir, target_suffix):
  for dirpath, dirnames, filenames in os.walk(root_dir):
    if 'openduck' in dirpath:
      for filename in filenames:
        full_file_path = os.path.join(dirpath, filename)
        if full_file_path.endswith(target_suffix):
          return os.path.abspath(full_file_path)
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
    replaceInFile(scriptPath, reps[0], reps[1])