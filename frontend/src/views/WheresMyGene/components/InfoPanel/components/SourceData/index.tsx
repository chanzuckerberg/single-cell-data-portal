import { List, ListItem } from "czifui";
import { useMemo } from "react";
import { Content, Header, ListSubheader, Wrapper } from "./style";

const MOCK_COLLECTION_URL = "https://czi.org";

const DATA_SETS = [
  {
    collection_label:
      "LungMAP — The genomic, epigenomic and biophysical cues controlling the emergence of the gas exchange niche in the lung",
    collection_url: MOCK_COLLECTION_URL,
    id: "Developmental single-cell atlas of the murine lung",
    label: "Developmental single-cell atlas of the murine lung",
  },
  {
    collection_label:
      "COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas",
    collection_url: MOCK_COLLECTION_URL,
    id: "Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients",
    label:
      "Large-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients",
  },
  {
    collection_label:
      "COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas",
    collection_url: MOCK_COLLECTION_URL,
    id: "Medium-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients",
    label:
      "Medium-scale single-cell analysis reveals critical immune characteristics of COVID-19 patients",
  },
  {
    collection_label:
      "A molecular cell atlas of the human lung from single cell RNA sequencing",
    collection_url: MOCK_COLLECTION_URL,
    id: "Krasnow Lab Human Lung Cell Atlas, 10X",
    label: "Krasnow Lab Human Lung Cell Atlas, 10X",
  },
  {
    collection_label: "Tabula Muris Senis",
    collection_url: MOCK_COLLECTION_URL,
    id: "All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - 10x",
    label:
      "All - A single-cell transcriptomic atlas characterizes ageing tissues in the mouse - 10x",
  },
  {
    collection_label: "Tabula Sapiens",
    collection_url: MOCK_COLLECTION_URL,
    id: "Tabula Sapiens - All Cells",
    label: "Tabula Sapiens - All Cells",
  },
  {
    collection_label:
      "LungMAP — Human data from a broad age healthy donor group",
    collection_url: MOCK_COLLECTION_URL,
    id: "Single-cell multiomic profiling of human lungs across age groups",
    label: "Single-cell multiomic profiling of human lungs across age groups",
  },
  {
    collection_label:
      // eslint-disable-next-line sonarjs/no-duplicate-string
      "Single-cell transcriptomics of human T cells reveals tissue and activation signatures in health and disease",
    collection_url: MOCK_COLLECTION_URL,
    id: "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Etiam tempus leo id quam commodo lacinia.",
    label:
      "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Etiam tempus leo id quam commodo lacinia.",
  },
  {
    collection_label:
      "Single-cell transcriptomics of human T cells reveals tissue and activation signatures in health and disease",
    collection_url: MOCK_COLLECTION_URL,
    id: "Lorem ipsum doladipiscing elit. Etiam tempus leo id quam commodo lacinia.",
    label:
      "Lorem ipsum doladipiscing elit. Etiam tempus leo id quam commodo lacinia.",
  },
  {
    collection_label:
      "Single-cell transcriptomics of human T cells reveals tissue and activation signatures in health and disease",
    collection_url: MOCK_COLLECTION_URL,
    id: "Lo id quam commodo lacinia.",
    label: "Lo id quam commodo lacinia.",
  },
];

interface Dataset {
  id: string;
  label: string;
  collection_label: string;
  collection_url: string;
}

interface Collection {
  name: string;
  url: string;
  datasets: { id: string; label: string }[];
}

interface Collections {
  [name: string]: Collection;
}

export default function SourceData(): JSX.Element {
  const collections: Collections = useMemo(() => {
    return getCollections(DATA_SETS);
  }, []);

  return (
    <Wrapper>
      <Header>Source Data</Header>
      <Content>
        <List ordered>
          {Object.values(collections).map(({ name, url, datasets }) => {
            return (
              <ListItem ordered key={name} fontSize="xxxs">
                <List
                  key={name}
                  subheader={
                    <a href={url} target="_blank" rel="noopener">
                      <ListSubheader>{name}</ListSubheader>
                    </a>
                  }
                >
                  {datasets.map(({ id, label }) => {
                    return (
                      <ListItem fontSize="xxxs" key={id}>
                        {label}
                      </ListItem>
                    );
                  })}
                </List>
              </ListItem>
            );
          })}
        </List>
      </Content>
    </Wrapper>
  );
}

function getCollections(datasets: Dataset[]): Collections {
  const collections: Collections = {};

  for (const dataset of datasets) {
    const { collection_label, collection_url, id, label } = dataset;

    if (!collections[collection_label]) {
      collections[collection_label] = {
        datasets: [],
        name: collection_label,
        url: collection_url,
      };
    }

    collections[collection_label].datasets.push({ id, label });
  }

  return collections;
}
