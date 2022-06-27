docker run 	--rm -it 	^
		-v "%CD%\work:/home/skimpy/work" 	^
		-v "%CD%/..:/skimpy"		^
		skimpy_docker %*
