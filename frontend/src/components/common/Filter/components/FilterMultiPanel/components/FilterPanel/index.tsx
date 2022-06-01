import { Divider, ListSubheader } from "@material-ui/core";
import { useEffect, useRef, useState } from "react";
import {
  CATEGORY_KEY,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
} from "src/components/common/Filter/common/entities";
import {
  Panel,
  useFilterPanelStyles,
} from "src/components/common/Filter/components/FilterMultiPanel/components/FilterPanel/style";
import FilterPanelList from "src/components/common/Filter/components/FilterMultiPanel/components/FilterPanelList";

const ADDITIONAL_PANEL_WIDTH = 8;

interface Props {
  categoryKey: CATEGORY_KEY;
  label?: string;
  onFilter: OnFilterFn;
  scrollable?: boolean;
  showPanelDivider: boolean;
  values: OntologyCategoryTreeNodeView[];
}

export default function FilterPanel({
  categoryKey,
  label,
  onFilter,
  scrollable = false,
  showPanelDivider,
  values,
}: Props): JSX.Element {
  const panelRef = useRef<HTMLDivElement>(null);
  const [panelWidth, setPanelWidth] = useState<number>(0);
  const classes = useFilterPanelStyles();

  // Calculate and set a min width on panel to prevent panel width resizing
  // (derived from a change in list item selected state font weight).
  useEffect(() => {
    if (panelRef.current) {
      setPanelWidth(
        panelRef.current?.children[0]?.clientWidth + ADDITIONAL_PANEL_WIDTH
      );
    }
  }, []);

  return (
    <>
      {showPanelDivider && (
        <Divider
          classes={{ root: classes.panelDivider }}
          orientation="vertical"
        />
      )}
      <Panel panelWidth={panelWidth} ref={panelRef} scrollable={scrollable}>
        <FilterPanelList
          categoryKey={categoryKey}
          onFilter={onFilter}
          PanelHeader={
            <ListSubheader
              classes={{ root: classes.panelHeading }}
              disableGutters
              disableSticky
              inset
            >
              {label}
            </ListSubheader>
          }
          values={values}
        />
      </Panel>
    </>
  );
}
