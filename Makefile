
.PHONY: environment remove-env install uninstall test # .PHONY is something we can add when our target dependencies are not files.

MODULE=mcsim
ENVIRONMENT=chem274A_lab3

environment: remove-env
	conda env create -f environment.yaml

remove-env:
	conda remove --name $(ENVIRONMENT) --all --yes

# Clean target to remove files except 'build_sys.sh' in 'octane' and 'diethylene_glycol' directories, but not in subdirectories
clean:
	@find octane diethylene_glycol -maxdepth 1 -type f ! -name 'build_sys.sh' -delete

