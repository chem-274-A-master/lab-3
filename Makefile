
.PHONY: environment remove-env clean build_sys # .PHONY is something we can add when our target dependencies are not files.

ENVIRONMENT=chem274A_lab3

environment: remove-env
	conda env create -f environment.yaml

remove-env:
	conda remove --name $(ENVIRONMENT) --all --yes

build_sys: clean
	@echo "Building system..."
	@cd build_box && bash build_sys.sh

# Clean target to remove files except 'build_sys.sh' in 'octane' and 'diethylene_glycol' directories, but not in subdirectories
clean:
	@find build_box -maxdepth 1 -type f ! -name 'build_sys.sh' -delete
	@find build_box/antechamber -type f ! -name '.gitkeep' -delete

