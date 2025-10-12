### Martini 3 development protein parameters ###
### Model #3                                 ###
### by L. Borges-Araujo, 2023-2024           ###

### CHANGELOG

#### [0.1.0.0.0] - Alpha release - February 2024 - LBA
##### Modified
 - Adjusted backbone masses to increase stability.
 - Adjusted stiff_fc from 100000 to 30000.
 - Adjusted alpha helix dihedral parameters. Now also using proper dihedrals like the remaining SS motifs.
 - Proline N partial charge moved to C. Mimics loss of amide hydrogen in proline. Reduces HBond estabilishment and reduces dipole.
 - Proline bonded parameters further refined from CHARMM, AMBER & PDB references.
 - Added specific dihedral parameters for proline (Coil & helix for now.)(WIP - Helix has only custom PSI.).
 - Added specific dihedral parameters for glycine.
 - ARG bonded parameters weakened.
##### Added
 - Introduced dFLEXIBLE flag to facilitate minimization.
 - Introduced new protein terminals. 
     - Neutral terminals -COOH -NH2.
     - Charged terminals.
 - Introduced amber forcefield mappings.
 - Introduced hydroxyproline. (WIP - No amber mappings yet. Bonded params same as PRO except for SC1 - BB bond which has been adjusted to account for hydroxyl.)
 - Introduced neutral GLY ASP LYS.
 - Added eh labels to TC5 beads in aromatic residues. Following directions from Ilias' and Liguo's work. Still under evaluation.
 
#### [0.0.3.9.4] - 2023-06-20 - LBA
##### Modified
 - Getting rid of -DSTABILIZE flag. All CA-SC1 bonds are now harmonic without need for flag.
 
#### [0.0.3.9.3] - 2023-06-12 - LBA
##### MAJOR
 - Improved COIL, BEND & TURN dihedrals. Should improve IDP & elastic network performance.
 - Fully revamped Martinize implementation. Fully fixes Martinize 2 issues.
##### Modified
 - Changed BB bead names to AA names. Makes most AA tools simply work.
 
#### [0.0.3.9.0] - 2023-05-03 - LBA
##### MAJOR
 - Substantially improved the bonded parameters of all sidechain residues.
##### Modified
 - Expanded LEU, ISO SC bond lengths by 15% to reduce stickyness.
 - Kept LEU and ISO as regular beads with r label.
 
#### [0.0.3.8.6_f] - 2023-04-26 - LBA
##### Reverted
 - Reverted on ISO, LEU size increase. Back to S SC beads.
##### Modified
 - Expanded LEU, ISO SC bond lengths by 15% to reduce stickyness.
 - Modified bonded parameters.
 
#### [0.0.3.8.6_e] - 2023-04-24 - LBA
##### Reverted
 - Reverted on ISO, LEU size increase. Back to S SC beads.
##### Modified
 - Expanded LEU, ISO SC bond lengths by 15% to reduce stickyness.
 
#### [0.0.3.8.6_d] - 2023-04-24 - LBA
##### Reverted
 - Reverted on r label on all C  beads. Kept only on SC.
##### Modified
 - Expanded LEU, ISO SC bond lengths by 15% to reduce stickyness. Kept R beads with r label.
 - Partial Charge set to 0.05.

#### [0.0.3.8.6_c] - 2023-04-24 - LBA
##### Reverted
 - Reverted on r label on all C  beads. Kept only on SC.
##### Modified
 - Expanded LEU, ISO SC bond lengths by 15% to reduce stickyness. Kept R beads with r label.
 - Partial Charge set to 0.0 (Effectively model #1).
 
#### [0.0.3.8.6_a] - 2023-04-20 - LBA
##### Reverted
 - Reverted on r label on all C  beads. Kept only on SC.
##### Modified
 - Expanded LEU, ISO SC bond lengths by 15% to reduce stickyness. Kept R beads with r label.
 
#### [0.0.3.8.6_b] - 2023-04-20 - LBA
##### Reverted
 - Reverted on r label on all C  beads. Kept only on SC.
 - Reverted on ISO, LEU size increase. Back to S SC beads.
 
#### [0.0.3.8.5] - 2023-04-11 - LBA
##### Modified
 - Adjusted SC1-CA -N  angle from 109.5 to 111.5 to match PDB distros.
 
#### [0.0.3.8.4] - 2023-04-11 - LBA
##### Added
 - Added r label to all C  beads to reduce protein-protein interactions.
 
#### [0.0.3.8.3] - 2023-04-05 - LBA
##### Added
 - First CHANGELOG for model #3.
 - Added r label to SC beads that had S->R size increase (LEU, ISO). For now, kept ASN without an r label. This can be revisited later.

