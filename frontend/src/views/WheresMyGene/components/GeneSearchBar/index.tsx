import { Intent } from "@blueprintjs/core";
import React, { useCallback, useContext, useMemo } from "react";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { usePrimaryFilterDimensions } from "src/common/queries/wheresMyGene";
import Toast from "src/views/Collection/components/Toast";
import { DispatchContext, StateContext } from "../../common/store";
import { selectGenes, selectTissues } from "../../common/store/actions";
import { Gene } from "../../common/types";
import GeneSets from "./components/Genesets";
import Organism from "./components/Organism";
import QuickSelect from "./components/QuickSelect";
import { ActionWrapper, Container } from "./style";

const GENESETS = [
  [
    "AGER",
    "FCN1",
    "CCL5",
    "PRF1",
    "BMX",
    "CCL23",
    "MS4A2",
    "RIMS2",
    "TP63",
    "FLT3",
    "PAX5",
  ],
  [
    "AGER",
    "FCN1",
    "CCL5",
    "PRF1",
    "BMX",
    "CCL23",
    "MS4A2",
    "RIMS2",
    "TP63",
    "FLT3",
    "PAX5",
    "MALAT1",
  ],
  [
    "LDB2",
    "VWF",
    "CA4",
    "PTPRB",
    "ADH1B",
    "GALNT18",
    "MAGI1",
    "KRT5",
    "TSPAN8",
    "ADRB1",
    "PLVAP",
    "PDGFRB",
    "MS4A1",
    "ACKR1",
    "RAMP3",
    "GNLY",
    "LTB",
    "GPR183",
    "PLEK",
    "TBXAS1",
    "AOAH",
    "ARHGAP15",
    "TPM2",
    "CALD1",
    "TACSTD2",
    "S100A8",
    "AIF1",
    "MS4A6A",
    "FGL2",
    "LYZ",
    "BTG1",
    "IL7R",
    "TAGLN",
    "ST6GALNAC5",
    "GPC5",
    "PDZRN3",
    "SFTA3_ENSG00000229415",
    "TP63",
    "LAMC3",
    "CSRP3-AS1",
    "LMNTD1",
    "GKN2",
    "PLA2G1B",
    "KRT23",
    "GABRP",
    "CFAP126",
    "LRRC10B",
    "FAM3D",
    "MUC4",
    "RTKN2",
    "SKAP1",
    "BLK",
    "SAMD3",
    "TPRG1",
    "DERL3",
    "MZB1",
    "CD68",
    "DEFB1",
    "HLA-DRB5",
    "CCL7",
    "HLA-DQA2",
    "STAC",
    "CP",
    "GRHL1",
    "MCEMP1",
    "TREM2",
    "RP11-1143G9.4",
    "S100A12",
    "CPA3",
    "TPSAB1",
    "TRBC2",
    "CD8A",
    "CD8B",
  ],
];

// DEBUG
// DEBUG
// DEBUG
// DEBUG

interface Tissue {
  name: string;
}

// END DEBUG

export default function GeneSearchBar(): JSX.Element {
  const dispatch = useContext(DispatchContext);
  const { selectedGenes, selectedTissues, selectedOrganismId } =
    useContext(StateContext);

  const { data } = usePrimaryFilterDimensions();

  const { genes: rawGenes, tissues } = data || {};

  const genes: Gene[] = useMemo(() => {
    if (!rawGenes) return [];

    return rawGenes[selectedOrganismId || ""] || [];
  }, [rawGenes, selectedOrganismId]);

  const genesByName = useMemo(() => {
    return genes.reduce((acc, gene) => {
      return acc.set(gene.name, gene);
    }, new Map<Gene["name"], Gene>());
  }, [genes]);

  const tissuesByName = useMemo(() => {
    const result = new Map<string, Tissue>();

    if (!tissues) return new Map<string, Tissue>();

    return tissues.reduce((acc, tissue) => {
      return acc.set(tissue.name, tissue);
    }, result);
  }, [tissues]);

  const selectedTissueOptions: Tissue[] = useMemo(() => {
    return selectedTissues.map((tissue: string) => {
      return tissuesByName.get(tissue) as Tissue;
    });
  }, [selectedTissues, tissuesByName]);

  const selectedGeneOptions: Gene[] = useMemo(() => {
    return selectedGenes.map((gene: string) => {
      return genesByName.get(gene) as Gene;
    });
  }, [selectedGenes, genesByName]);

  const handleGeneNotFound = useCallback((geneName: string): void => {
    Toast.show({
      intent: Intent.DANGER,
      message: `Gene not found: ${geneName}`,
    });
  }, []);

  return (
    <Container>
      {/* DEMO ONLY WILL BE DELETED BEFORE MVP */}
      {/* DEMO ONLY WILL BE DELETED BEFORE MVP */}
      {/* DEMO ONLY WILL BE DELETED BEFORE MVP */}
      <GeneSets onSelect={handleGenesetsSelect} />

      <br />
      <br />
      <ActionWrapper>
        <Organism />

        <QuickSelect
          items={tissues || EMPTY_ARRAY}
          itemsByName={tissuesByName}
          multiple
          selected={selectedTissueOptions}
          setSelected={handleSelectTissues}
          label="Add Tissue"
          dataTestId="add-tissue"
        />

        <QuickSelect
          items={genes}
          itemsByName={genesByName}
          selected={selectedGeneOptions}
          multiple
          setSelected={handleSelectGenes}
          onItemNotFound={handleGeneNotFound}
          label="Add Gene"
          dataTestId="add-gene"
        />
      </ActionWrapper>
    </Container>
  );

  function handleGenesetsSelect(genesetIndex: number) {
    const geneset = GENESETS[genesetIndex];

    handleSelectGenes(
      genes
        .filter((gene) => geneset.includes(gene.name))
        .sort(
          (a, b) =>
            geneset.findIndex((gene) => gene === a.name) -
            geneset.findIndex((gene) => gene === b.name)
        )
    );
  }

  function handleSelectTissues(tissues: Tissue[]) {
    if (!dispatch) return;

    dispatch(selectTissues(tissues.map((tissue) => tissue.name)));
  }

  function handleSelectGenes(genes: Gene[]) {
    if (!dispatch) return;

    dispatch(selectGenes(genes.map((gene) => gene.name)));
  }
}
