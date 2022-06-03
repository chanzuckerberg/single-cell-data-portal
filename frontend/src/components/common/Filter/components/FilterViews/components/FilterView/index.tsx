import { Divider, ListSubheader } from "@material-ui/core";
import { useEffect, useRef, useState } from "react";
import {
  CATEGORY_KEY,
  OnFilterFn,
  OntologyCategoryTreeNodeView,
} from "src/components/common/Filter/common/entities";
import {
  useFilterViewStyles,
  ViewPanel,
} from "src/components/common/Filter/components/FilterViews/components/FilterView/style";
import FilterViewList from "src/components/common/Filter/components/FilterViews/components/FilterViewList";

const ADDITIONAL_PANEL_WIDTH = 8;

interface Props {
  categoryKey: CATEGORY_KEY;
  label?: string;
  onFilter: OnFilterFn;
  scrollable?: boolean;
  showViewDivider: boolean;
  values: OntologyCategoryTreeNodeView[];
}

export default function FilterView({
  categoryKey,
  label,
  onFilter,
  scrollable = false,
  showViewDivider,
  values,
}: Props): JSX.Element {
  const panelRef = useRef<HTMLDivElement>(null);
  const [panelWidth, setPanelWidth] = useState<number>(0);
  const classes = useFilterViewStyles();

  // Calculate and set a min width on menu panel to prevent width resizing
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
      {showViewDivider && (
        <Divider
          classes={{ root: classes.viewDivider }}
          orientation="vertical"
        />
      )}
      <ViewPanel panelWidth={panelWidth} ref={panelRef} scrollable={scrollable}>
        <FilterViewList
          categoryKey={categoryKey}
          onFilter={onFilter}
          values={values}
          ViewHeader={
            <ListSubheader
              classes={{ root: classes.viewHeading }}
              disableGutters
              disableSticky
              inset
            >
              {label}
            </ListSubheader>
          }
        />
      </ViewPanel>
    </>
  );
}
