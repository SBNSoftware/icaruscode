ICARUS CVN Mapper and Evaluator
Author: Felipe A. Wieler
Date: November 12, 2025
Email: f3aw@proton.me

------------------------------------------------------------
Overview
------------------------------------------------------------
The ICARUSCVNMapper and ICARUSCVNEvaluator modules provide the interface between reconstructed slices and the Convolutional Visual Network (CVN) models used in ICARUS.

- ICARUSCVNMapper.fcl: Configuration file for detector mapping and label management, including Grid Configuration.
- ICARUSZlibMaker.fcl: Converts the mapper output into compressed pixel map archives.
- ICARUSCVNEvaluator.fcl: Runs the trained CVN TensorFlow models to evaluate each slice and produce classification scores.

These tools together enable pixel map creation, compression, and evaluation for neutrino flavors classification in ICARUS data.

------------------------------------------------------------
Prerequisites
------------------------------------------------------------
Before running, ensure the following:
- LArSoft and ICARUS code are correctly set up in your environment.
- TensorFlow model files are accessible, those are defined on the ICARUSCVNMapper.fcl.
- Input files (stage1) have been produced by standard ICARUS reconstruction chains.
- standard_cvnzlibmaker is properly configured, including the OutputDir path also in ICARUSCVNMapper.fcl.

------------------------------------------------------------
Pixel Map Generation
------------------------------------------------------------
To create pixel maps from a stage1 file:
lar -c run_ICARUSCVNMapper.fcl -s <stage1_file.root>

This produces a ROOT output with slice and pixel map information.

Next, run the zlib maker to generate the compressed maps:
lar -c run_ICARUSZlibMaker.fcl -s <mapper_output.root>

This step creates:
- .gz files → contain pixel map arrays
- .info files → include truth and slice metadata from Pandora

These files are stored under the directory specified by OutputDir in standard_cvnzlibmaker on ICARUSCVNMapper.fcl.

------------------------------------------------------------
Running the Evaluator
------------------------------------------------------------
To evaluate slices using the trained CVN models:
lar -c run_ICARUSCVNEvaluator.fcl -s <stage1_file.root>

This produces a new ROOT file containing the CryoE and CryoW classification scores for each slice.

------------------------------------------------------------
Output Summary
------------------------------------------------------------
File Type | Description
----------|------------------------------------------------
.gz       | Compressed pixel map arrays
.info     | Truth and metadata from Pandora per slice
.root     | Evaluation results including CryoE and CryoW scores

------------------------------------------------------------
Example Workflow
------------------------------------------------------------

Pixel Map Creation:
lar -c run_ICARUSCVNMapper.fcl     -s reco_stage1.root
lar -c run_ICARUSZlibMaker.fcl     -s output_mapper.root
------------------------------------------------------------

Root Evaluation:
lar -c run_ICARUSCVNEvaluator.fcl  -s reco_stage1.root

------------------------------------------------------------
Notes
------------------------------------------------------------
- Use consistent stage1 files across all steps to maintain slice–truth consistency.
- Evaluation requires access to trained CVN TensorFlow models.
- The modules support both Cryostat East (CryoE) and Cryostat West (CryoW) toghether or individualy by changin the corresponding FHiCL file.
