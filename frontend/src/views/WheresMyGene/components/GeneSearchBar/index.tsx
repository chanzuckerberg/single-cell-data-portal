import { Intent } from "@blueprintjs/core";
import { ButtonBase, Popper, Theme } from "@material-ui/core";
import { AutocompleteCloseReason } from "@material-ui/lab";
import { makeStyles } from "@material-ui/styles";
import {
  DefaultMenuSelectOption,
  getColors,
  getCorners,
  getShadows,
  MenuSelect,
} from "czifui";
import pull from "lodash/pull";
import uniq from "lodash/uniq";
import React, {
  createContext,
  useEffect,
  useMemo,
  useRef,
  useState,
} from "react";
import { FixedSizeList, ListChildComponentProps } from "react-window";
import { API } from "src/common/API";
import { EMPTY_ARRAY } from "src/common/constants/utils";
import { DEFAULT_FETCH_OPTIONS } from "src/common/queries/common";
import { API_URL } from "src/configs/configs";
import Toast from "src/views/Collection/components/Toast";
import { Gene } from "../../common/types";
import GeneSets from "./components/Genesets";
import { Container } from "./style";

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

const LISTBOX_ITEM_HEIGHT_PX = 32;
const LISTBOX_HEIGHT_PX = 152;

interface Props {
  onGenesChange: (selectedGenes: Gene[]) => void;
}

const ListBoxContext = createContext({});

const OuterElementType = React.forwardRef<HTMLDivElement>(
  function OuterElementType(props, ref) {
    const outerProps = React.useContext(ListBoxContext);
    return <div ref={ref} {...props} {...outerProps} />;
  }
);

function rowRender(props: ListChildComponentProps) {
  const { data, index, style } = props;
  return <div style={style}>{data[index]}</div>;
}

const ListboxComponent = React.forwardRef<HTMLDivElement>(
  function ListboxComponent(props, ref) {
    const { children, ...other } = props;

    const itemData = React.Children.toArray(children);
    const itemCount = itemData.length;

    return (
      <div ref={ref}>
        <ListBoxContext.Provider value={other}>
          <FixedSizeList
            height={LISTBOX_HEIGHT_PX}
            itemCount={itemCount}
            outerElementType={OuterElementType}
            itemSize={LISTBOX_ITEM_HEIGHT_PX}
            width="100%"
            overscanCount={10}
            itemData={itemData}
          >
            {rowRender}
          </FixedSizeList>
        </ListBoxContext.Provider>
      </div>
    );
  }
);

export default function GeneSearchBar({ onGenesChange }: Props): JSX.Element {
  const [selectedGenes, setSelectedGenes] = useState<Gene[]>(EMPTY_ARRAY);
  const [genes, setGenes] = useState<Gene[]>(EMPTY_ARRAY);
  const [open, setOpen] = useState(false);
  const [pendingPaste, setPendingPaste] = useState(false);
  const [input, setInput] = useState("");

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

  const handleClose = (
    _: React.ChangeEvent<Record<string, never>>,
    reason: AutocompleteCloseReason
  ) => {
    if (reason === "toggleInput") {
      return;
    }
    setOpen(false);
  };
  const handleChange = (
    _: React.ChangeEvent<Record<string, never>>,
    newValue: DefaultMenuSelectOption[] | null
  ) => {
    return setSelectedGenes(newValue as Gene[]);
  };
  const handleClick = () => {
    setOpen(true);
  };

  const handlePaste = () => {
    setPendingPaste(true);
  };

  const genesByName = useMemo(() => {
    return genes.reduce((acc, gene) => {
      return acc.set(gene.name, gene);
    }, new Map<Gene["name"], Gene>());
  }, [genes]);

  const handleEnter = (event: React.KeyboardEvent<HTMLInputElement>) => {
    if (event.key === "Enter" && pendingPaste) {
      event.preventDefault();
      const newSelectedGenes = [...selectedGenes];
      const pastedGenes = pull(uniq(input.split(/[ ,]+/)), "");
      pastedGenes.map((gene) => {
        const newGene = genesByName.get(gene);
        if (!newGene) {
          Toast.show({
            intent: Intent.DANGER,
            message: `Gene not found: ${gene}`,
          });
        } else if (!newSelectedGenes.includes(newGene))
          newSelectedGenes.push(newGene);
      });
      setPendingPaste(false);
      setOpen(false);
      return setSelectedGenes(newSelectedGenes);
    }
  };

  const useStyles = makeStyles((theme: Theme) => {
    const colors = getColors({ theme });
    const shadows = getShadows({ theme });
    const corners = getCorners({ theme });
    return {
      paper: {
        boxShadow: "none",
        margin: 0,
      },
      popper: {
        backgroundColor: "white",
        border: `1px solid ${colors?.gray[100]}`,
        borderRadius: corners?.m,
        boxShadow: shadows?.m,
        color: "#586069",
        fontSize: 13,
        width: 377,
        zIndex: 3, // The x axis wrapper is set at 2
      },
      popperDisablePortal: {
        position: "relative",
        width: "100% !important",
      },
    };
  });

  const classes = useStyles();

  const ref = useRef(null);

  return (
    <Container>
      <GeneSets onSelect={handleGenesetsSelect} />

      <br />
      <br />
      <ButtonBase disableRipple onClick={handleClick} ref={ref}>
        <span>Add Genes</span>
      </ButtonBase>

      <Popper open={open} className={classes.popper} anchorEl={ref.current}>
        <MenuSelect
          open
          search
          onClose={handleClose}
          multiple
          classes={{
            paper: classes.paper,
            popperDisablePortal: classes.popperDisablePortal,
          }}
          value={selectedGenes}
          onChange={handleChange}
          disableCloseOnSelect
          disableListWrap
          onKeyDownCapture={handleEnter}
          options={genes}
          ListboxComponent={
            ListboxComponent as React.ComponentType<
              React.HTMLAttributes<HTMLElement>
            >
          }
          renderOption={(option) => option.name}
          onPaste={handlePaste}
          // eslint-disable-next-line @typescript-eslint/ban-ts-comment -- TODO(cc) revisit
          // @ts-ignore -- TODO(cc) revisit
          InputBaseProps={{
            onChange: (
              event: React.ChangeEvent<HTMLTextAreaElement | HTMLInputElement>
            ) => {
              setInput(event.target.value);
            },
            placeholder: "Search or paste comma separated gene names",
          }}
        />
      </Popper>
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
}
