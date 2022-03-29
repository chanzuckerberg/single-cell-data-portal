import { List, ListItem } from "czifui";
import { useMemo } from "react";
import { ROUTES } from "src/common/constants/routes";
import {
  FilterDimensions,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { Content, Header, ListSubheader, Wrapper } from "./style";

interface Collection {
  name: string;
  url: string;
  datasets: { id: string; label: string }[];
}

interface Collections {
  [name: string]: Collection;
}

export default function SourceData(): JSX.Element {
  const { data: filterDimensions } = useFilterDimensions();
  const { datasets = [] } = filterDimensions;

  const collections: Collections = useMemo(() => {
    return getCollections(datasets);
  }, [datasets]);

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

function getCollections(datasets: FilterDimensions["datasets"]): Collections {
  const collections: Collections = {};

  for (const dataset of datasets) {
    const { collection_label, collection_id, id, label } = dataset;

    if (!collections[collection_label]) {
      collections[collection_label] = {
        datasets: [],
        name: collection_label,
        url: ROUTES.COLLECTION.replace(":id", collection_id),
      };
    }

    collections[collection_label].datasets.push({ id, label });
  }

  return collections;
}
