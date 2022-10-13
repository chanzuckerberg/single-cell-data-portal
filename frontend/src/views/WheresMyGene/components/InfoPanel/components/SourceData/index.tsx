import { List, ListItem } from "czifui";
import { useMemo } from "react";
import {
  aggregateCollectionsFromDatasets,
  useFilterDimensions,
} from "src/common/queries/wheresMyGene";
import { Content, ListSubheader, Wrapper } from "./style";

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
    return aggregateCollectionsFromDatasets(datasets);
  }, [datasets]);

  return (
    <Wrapper>
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
