import { Menu, MenuItem } from "@blueprintjs/core";
import {
  IItemListRendererProps,
  IItemRendererProps,
  MultiSelect,
} from "@blueprintjs/select";
import { forwardRef, useEffect, useState } from "react";
import { FixedSizeList } from "react-window";
import { API } from "src/common/API";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { API_URL } from "src/configs/configs";
import { Gene } from "../../common/types";
import GeneSets from "./components/Genesets";
import { Container, MultiSelectWrapper } from "./style";

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

interface Props {
  onGenesChange: (selectedGenes: Gene[]) => void;
}

interface ExtendedItemRendererProps extends IItemRendererProps {
  item: Gene;
  isSelected: boolean;
}

export default function GeneSearchBar({ onGenesChange }: Props): JSX.Element {
  const [selectedGenes, setSelectedGenes] = useState<Gene[]>(EMPTY_ARRAY);
  const [genes, setGenes] = useState<Gene[]>(EMPTY_ARRAY);
  // DEBUG
  // DEBUG
  // DEBUG
  // (thuang): Default to 40 genes for now
  const [input, setInput] = useState("40");

  useEffect(() => {
    fetchGenes();

    async function fetchGenes(): Promise<void> {
      const response = await fetch(
        API_URL + API.WMG_GENES,
        DEFAULT_FETCH_OPTIONS
      );

      // DEBUG
      // DEBUG
      // DEBUG
      // (thuang): Local test data
      // const response = await fetch(
      //   "https://wmg-prototype-data-dev-public.s3.amazonaws.com/lung-tissue-10x-human/lung_tissue_genes.json"
      // );

      const allGenes = await response.json();

      setGenes(allGenes);
    }
  }, []);

  useEffect(() => {
    onGenesChange(selectedGenes);
  }, [onGenesChange, selectedGenes]);

  useEffect(() => {
    if (Number.isNaN(Number(input))) return;

    setSelectedGenes(genes.slice(0, Number(input)));
  }, [input, genes]);

  return (
    <Container>
      <GeneSets onSelect={handleGenesetsSelect} />

      <br />
      <br />

      <label htmlFor="first-n-genes">Select first N genes</label>
      <br />
      <input id="first-n-genes" onChange={handleInputChange} value={input} />
      <MultiSelectWrapper>
        <MultiSelect
          itemPredicate={itemPredicate}
          onItemSelect={handleItemSelect}
          onRemove={handleItemRemove}
          items={genes}
          itemRenderer={renderItem}
          tagRenderer={TagRenderer}
          itemsEqual={areGenesEqual}
          selectedItems={selectedGenes}
          itemListRenderer={itemListRenderer}
        />
      </MultiSelectWrapper>
    </Container>
  );

  function handleGenesetsSelect(genesetIndex: number) {
    const geneset = GENESETS[genesetIndex];

    setSelectedGenes(
      genes
        .filter((gene) => geneset.includes(gene.name))
        .sort(
          (a, b) =>
            geneset.findIndex((gene) => gene === a.name) -
            geneset.findIndex((gene) => gene === b.name)
        )
    );
  }

  function handleInputChange(event: React.ChangeEvent<HTMLInputElement>): void {
    setInput(event.target.value);
  }

  function itemListRenderer(listProps: IItemListRendererProps<Gene>) {
    const {
      filteredItems,
      renderItem: propRenderItem,
      itemsParentRef,
    } = listProps;

    return (
      <FixedSizeList
        // eslint-disable-next-line react/display-name
        innerElementType={forwardRef((props, ref) => {
          return <Menu ulRef={itemsParentRef} ref={ref} {...props} />;
        })}
        height={300}
        overscanCount={5}
        width="100%"
        itemCount={filteredItems.length}
        itemSize={24}
      >
        {(props) => {
          const { index, style } = props;

          return propRenderItem({ ...filteredItems[index], style }, index);
        }}
      </FixedSizeList>
    );
  }

  function handleItemSelect(gene: Gene) {
    if (isGeneSelected(gene)) {
      handleItemRemove(gene);
    } else {
      setSelectedGenes((prevSelectedGenes) => [...prevSelectedGenes, gene]);
    }
  }

  function handleItemRemove(gene: Gene) {
    setSelectedGenes(
      selectedGenes.filter((selectedGene) => selectedGene.id !== gene.id)
    );
  }

  function renderItem(gene: Gene, itemRendererProps: IItemRendererProps) {
    return ItemRenderer({
      isSelected: isGeneSelected(gene),
      item: gene,
      ...itemRendererProps,
    });
  }

  function isGeneSelected(gene: Gene): boolean {
    return Boolean(
      selectedGenes.find((selectedGene) => selectedGene.id === gene.id)
    );
  }
}

function ItemRenderer({
  item,
  handleClick,
  query,
  isSelected,
}: ExtendedItemRendererProps): JSX.Element | null {
  const { name, style } = item;

  return (
    <MenuItem
      active={isSelected}
      key={name}
      onClick={handleClick}
      text={highlightText(name, query)}
      style={{ ...style }}
    />
  );
}

function itemPredicate(query: string, item: Gene) {
  return item.name.toLowerCase().indexOf(query.toLowerCase()) >= 0;
}

function highlightText(text: string, query: string) {
  let lastIndex = 0;
  const words = query
    .split(/\s+/)
    .filter((word) => word.length > 0)
    .map(escapeRegExpChars);
  if (words.length === 0) {
    return [text];
  }
  const regexp = new RegExp(words.join("|"), "gi");
  const tokens: React.ReactNode[] = [];
  // eslint-disable-next-line no-constant-condition -- expected use
  while (true) {
    const match = regexp.exec(text);
    if (!match) {
      break;
    }
    const wordLength = match[0].length;
    const before = text.slice(lastIndex, regexp.lastIndex - wordLength);
    if (before.length > 0) {
      tokens.push(before);
    }
    lastIndex = regexp.lastIndex;
    tokens.push(<strong key={lastIndex}>{match[0]}</strong>);
  }
  const rest = text.slice(lastIndex);
  if (rest.length > 0) {
    tokens.push(rest);
  }
  return tokens;
}

function escapeRegExpChars(text: string) {
  return text.replace(/([.*+?^=!:${}()|[\]/\\])/g, "\\$1");
}

function TagRenderer({ name }: Gene) {
  return name;
}

function areGenesEqual(geneA: Gene, geneB: Gene) {
  // Compare only the names (ignoring case) just for simplicity.
  return geneA.id.toLowerCase() === geneB.id.toLowerCase();
}
