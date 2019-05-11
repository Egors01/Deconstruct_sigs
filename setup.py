
from setuptools import find_packages, Distribution
from setuptools import setup

def get_module_name(module):
    return module.__name__.split('.')[-1]


class BinaryDistribution(Distribution):
    def is_pure(self):
        return False
#
# def package_files(directory):
#     import os
#     paths = []
#     for (path, directories, filenames) in os.walk(directory):
#         for filename in filenames:
#             if os.path.splitext(filename)[-1] != '.py':
#                 paths.append(os.path.join(path, filename))
#     return paths
# extra_files = package_files('deconstructsigs')
setup(
    name='deconstructsigs',
    version='0.3',
    packages=find_packages(),
    license='https://github.com/Egors01/deconstruct_sigs',
    author='egors',
    author_email='none',
    distclass=BinaryDistribution,
    description=''
)

#
# setup(name='deconstructsigs',
#       description='Implementation of DeconstructSigs_old algorithm for deducing cancer genome mutational signatures',
#       author='Eric Kofman mod Egors',
#       author_email='',
#       version='0.1',
#       py_modules=['deconstructsigs'],
#       url='https://github.com/Egors01/deconstruct_sigs',)
#


# def package_files(directory):
#     import os
#     paths = []
#     for (path, directories, filenames) in os.walk(directory):
#         for filename in filenames:
#             if os.path.splitext(filename)[-1] != '.py':
#                 paths.append(os.path.join(path, filename))
#     return paths
# extra_files = package_files('primatestools')
# # # primary_extra_files = [extra.replace('biota/', '') for extra in primary_extra_files]
# #
# # package_data.extend(extra_files)
# # print(package_data)
# #print (find_packages())
# # setup(
# #     name='primatestools',
# #     version='0.5',
# #     packages=find_packages(),
# #     url='https://github.com/Egors01/Mut_spectra_primates',
# #     license='',
# #     author='egors',
# #     author_email='none',
# #     distclass=BinaryDistribution,
# #     description=''
# # )