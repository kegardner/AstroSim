import os
import sys
import scopesim as sim

local_package_folder = "../inst_pkgs"
sim.rc.__config__["!SIM.file.local_packages_path"] = local_package_folder

packages = ["MICADO","Armazones", "ELT", "MORFEO"]

toDownload = list(filter(lambda p : not os.path.isdir(f"{local_package_folder}/{p}") or sys.argv.count("--update-ss"), packages))
#guard with an if as scopesim fires a few http calls even if the input is empty; plus it throws an exception when the list if empty
if any(toDownload) :
    sim.download_packages(toDownload) 

