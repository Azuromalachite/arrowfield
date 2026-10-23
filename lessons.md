# Organic-chemistry companion: lessons

Six short lessons in the order a first organic course (McMurry) meets the
ideas: polar bonds, carbonyls, acids and bases, esters, alkenes and
stereochemistry, resonance. Each uses molecules from `molecules.yaml`, runs
live in seconds, and follows one shape:

> **Predict** (pick an answer before anything is shown) → **Reveal** (surface
> + 2D structure + site table) → **Explain** (with ordinary organic reasoning)
> → **Control** (a molecule that isolates one change) → **Limits** (what this
> picture cannot tell you)

The answers below are what the chemistry leads you to expect. They are not
measured espmap numbers. Each lesson's reveal shows its own values, and where
a number is needed later it comes from a computed record.

Each lesson names the McMurry sections it accompanies as `McMurry §x.y`, to
be filled in for the class edition. Read the textbook there; the lesson adds
only the computed picture. Chemistry claims cite `docs/references.bib` keys,
or are marked `TODO: cite`.

---

## 1. Polar covalent bonds

**Molecules:** ethane (control) → ethanol → dimethyl ether · **Read:** McMurry §__ (polar covalent bonds, electronegativity) · `politzer2002`

**Predict**
- Where is ethanol most electron-rich? ( ) the O ( ) the C bonded to O ( ) the CH₃ ( ) not sure
- Which hydrogen is the most electron-poor? ( ) the O–H ( ) a C–H ( ) all the same

**Reveal:** sites `hydroxyl` (O + its H) and `C-O carbon` (C2 + H).

**Explain**
- Electronegativity: O pulls density from C and H, so the C–O and O–H bonds are
  polar. The lone pairs make the region around O negative.
- The O–H hydrogen is strongly positive: this is the hydrogen-bond donor.
  The lone pairs on O are the acceptor.

**Control:** dimethyl ether has the same formula, C₂H₆O. Its O is still
negative, but there is no positive O–H: an acceptor without a donor. Ethane
shows what "no polar bond" looks like on the same scale.

**Limits:** a positive region on H says where a negative partner is
attracted. It does not give the pKa: acidity depends on the conjugate base
(lesson 3).

---

## 2. The carbonyl group: red is not "reactive"

**Molecules:** ethanol → acetaldehyde → acetone · **Read:** McMurry §__ (the carbonyl group; nucleophilic addition)

**Predict**
- Where is acetaldehyde most electron-rich? ( ) carbonyl O ( ) carbonyl C ( ) the CH₃
- A nucleophile attacks a carbonyl at ( ) O ( ) C, so which matters for that, the red region or the C next to it?

**Reveal:** sites `carbonyl O` and `carbonyl C-H` (aldehyde) or `carbonyl C` (ketone).

**Explain**
- C=O is polarised, Cδ+–Oδ−. The negative region sits on O, and the carbonyl C
  is comparatively electron-poor.
- A nucleophile attacks the **carbonyl C**, not the red O. The surface shows
  the polarity that makes the C electrophilic; it does not point at the
  reacting atom. This is the first place "red = reactive" fails.
- Acetone's two methyl groups donate density to the C=O carbon. Compare its
  carbonyl C with acetaldehyde's on one scale (same `render_key`).

**Control:** ethanol, the reduced form: C–O single bond, no π bond, no
electrophilic carbonyl carbon.

**Limits:** which carbonyl reacts faster is about barriers and sterics, not
just surface polarity.

---

## 3. Acids and their conjugate bases

**Molecules:** acetic acid → acetate (**PCM water for both**) · **Read:** McMurry §__ (acids and bases; resonance)

**Geometry:** acetate needs a DFT geometry optimised in PCM water before this
lesson runs (manifest `TODO`). SMILES + MMFF can move a small anion by up to
6.5 kcal/mol (H8).

**Predict**
- In acetate, are the two oxygens ( ) different, one C=O and one C–O⁻, as drawn ( ) the same
- The whole acetate surface is ( ) mostly negative ( ) half positive, half negative

**Reveal:** the colour scale switches to **sequential**: an ion's surface
never changes sign, and espmap refuses to centre it on zero.

**Explain**
- Resonance: the negative charge is shared over both oxygens. The drawn
  structure is one resonance form, and the surface shows the average. The
  two O sites come out the same (espmap's symmetry test checks this on
  nitrite to 0.05 kcal/mol).
- Delocalisation stabilises the conjugate base, which is why carboxylic acids
  are more acidic than alcohols.

**Control:** acetic acid in the same solvent and on a shared scale. The
O–H is strongly positive, and the C=O and O–H oxygens are not equivalent.

**Limits:**
- An ion and a neutral molecule cannot share one diverging scale.
- A gas-phase anion is not a stand-in for the solvated one. Here, gas-phase
  small anions are actually unbound.

---

## 4. Esters: from the textbook to the pheromone

**Molecules:** ethanol → acetic acid → ethyl acetate · **Read:** McMurry §__ (carboxylic acid derivatives; esters)

**Predict**
- Which ester oxygen is more electron-rich? ( ) the carbonyl O ( ) the ether-like O ( ) the same

**Reveal:** sites `carbonyl O` and `ester O`.

**Explain**
- The ether-type O donates a lone pair into the C=O by resonance, so the
  carbonyl O carries most of the negative potential.
- An ester has no O–H donor, so the acceptor-only pattern of lesson 1's ether
  returns.
- This is the acetate head group of the Sesiidae pheromones: the same local
  chemistry sits on the end of an 18-carbon chain.

**Control:** acetic acid, which is the same carbonyl but with a donor O–H.

**Limits:** the chain contributes almost nothing to the head group's local
polarity. Its influence on how the molecule is recognised is shape, which
comes in lesson 5 and the homework.

---

## 5. Alkenes and E/Z

**Molecules:** ethene → (E)-but-2-ene → (Z)-but-2-ene · **Read:** McMurry §__ (alkene structure; E,Z designation)

**Predict**
- Is an alkene's π bond ( ) electron-rich above and below the plane ( ) electron-poor ( ) neutral
- Do E- and Z-but-2-ene differ ( ) mainly in charge distribution ( ) mainly in shape

**Reveal:** both butenes on one shared scale, same camera.

**Explain**
- The π electrons sit above and below the C=C plane: a mildly negative
  region, which is why alkenes react with electrophiles.
- E and Z have the same atoms and bonds, and differ in where the methyls are.
  Expect the local π region to look alike and the overall shape to differ.
  The homework asks the same question of the C18 pheromones, where it is
  harder.

**Control:** ethene, the C=C with no substituents.

**Limits:** how much E/Z changes a number has to be compared with how much
the same molecule's conformers change it. That is the difference card, and
the homework's whole point.

---

## 6. Resonance and two reactive sites (an enolate)

**Molecules:** acetaldehyde → acetaldehyde enolate (**PCM water**) · **Read:** McMurry §__ (enolate ions; alkylation of enolates)

**Geometry:** espmap's `bench/fixtures/enolate.sdf` (B3LYP/def2-TZVP, PCM
water), never SMILES + MMFF.

**Predict**
- The enolate's negative charge is on ( ) O only ( ) C only ( ) shared by O and the CH₂ carbon
- An alkyl halide (SN2) attacks ( ) O ( ) C ( ) it depends. If it depends, on what?

**Reveal:** sites O and the CH₂ carbon.

**Explain**
- Two resonance forms put the negative charge on O or on the CH₂ carbon: an
  ambident nucleophile.
- On the surface O is clearly the most negative site: the charge-controlled
  site. The calibration pilot found the same for all three ambident anions
  it tried (enolate O, thiocyanate N, nitrite O).
- **Where the anion reacts is not one fact about the anion.** It depends on
  the electrophile, the counter-ion, the solvent, and whether the fastest
  product or the most stable one wins:
  - **The enolate.** A free enolate is alkylated faster at O (the kinetic
    product); the C-alkylated product is the more stable one (the
    thermodynamic product) `mayr2011` (§3.4.2). The C-alkylation that
    synthesis relies on comes from metal-enolate aggregates, the
    counter-ion and the solvent, none of which is in a picture of the
    free ion.
  - **Nitrite.** Methylation gives O:N from 30:70 to 92:8 as the solvent
    and the methylating agent change `tishkov2005` (Table 2). Sodium nitrite
    with primary and secondary alkyl bromides in DMF gives mostly the nitro
    compound (55–65 %), with 30–33 % of the nitrite ester `kornblum1955`.
  - **Thiocyanate** is the clean case: SN2 alkylation goes to S, favoured
    over N by 10²–10³ `loos2003`. The most negative site on its surface is
    N, so here the surface and the reaction really do disagree.
- In the pilot, ALIE (a descriptor of how easily electrons are removed)
  marked thiocyanate's S, and the enolate's CH₂ carbon, the site of the more
  stable product. For nitrite neither ALIE nor the orbital descriptors
  separated N from O clearly, which fits an anion whose selectivity changes
  with the conditions.
- The lesson: the electrostatic picture shows charge control. Where a
  reaction happens can be decided by something else: orbitals,
  transition-state energies, and the conditions, which a calculation on the
  isolated anion does not contain. A textbook's "alkylates at C" summarises
  particular conditions; find out which.

**Control:** acetaldehyde, the neutral carbonyl it comes from.

**Extension (solvent):** phenoxide alkylates at O in DMF, THF and alcohols,
but substantially at C in water or phenol, attributed to hydrogen bonding at
O `kornblum1963`. The same anion gives a different answer, decided by how
the solvent holds its oxygen.

**Limits:** this is where "ESP shows the polarity" and "ESP predicts the
reaction" come apart, and where later tools (ALIE, frontier orbitals) earn
their place. Even they describe the isolated anion; the counter-ion and the
solvent can decide the outcome.

Sources for this lesson's outcomes (keys from espmap's `docs/references.bib`,
added in espmap #66):
[`mayr2011`](https://doi.org/10.1002/anie.201007100),
[`tishkov2005`](https://doi.org/10.1002/anie.200501274),
[`kornblum1955`](https://doi.org/10.1021/ja01628a064),
[`loos2003`](https://doi.org/10.1021/ja037317u),
[`kornblum1963`](https://doi.org/10.1021/ja00891a024).

---

## Pacing

- One lesson per sitting.
- Lessons 1–4 need only bond polarity and resonance.
- Lesson 5 needs E/Z naming.
- Lesson 6 is best after nucleophilic substitution has been met.
- The homework assumes lessons 1, 2, 4 and 5.
