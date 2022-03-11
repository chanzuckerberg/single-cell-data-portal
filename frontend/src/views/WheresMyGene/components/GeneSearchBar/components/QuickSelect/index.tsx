import { Theme } from "@emotion/react";
import { makeStyles, Popper } from "@material-ui/core";
import { AutocompleteCloseReason } from "@material-ui/lab";
import {
  DefaultMenuSelectOption,
  getColors,
  getCorners,
  getShadows,
  Icon,
  MenuSelect,
} from "czifui";
import { pull, uniq } from "lodash";
import React, { createContext, useRef, useState } from "react";
import { FixedSizeList, ListChildComponentProps } from "react-window";
import { noop } from "src/common/constants/utils";
import { Label } from "../../style";
import { ButtonWrapper, StyledIconButton } from "./style";

const LISTBOX_ITEM_HEIGHT_PX = 32;
const LISTBOX_HEIGHT_PX = 152;

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

// (thuang): Value's type is based on generic type placeholder (T) and Multiple
// type. If Multiple is true, Value's type is T[].
// Otherwise, Value's type is T.
// Conditional Type
// https://www.typescriptlang.org/docs/handbook/2/conditional-types.html
export type Value<T, Multiple> = Multiple extends undefined | false ? T : T[];

interface Props<T, Multiple> {
  items: T[];
  multiple?: Multiple;
  setSelected: (selected: Value<T, Multiple>) => void;
  selected: Value<T, Multiple>;
  itemsByName: Map<string, T>;
  onItemNotFound?: (item: string) => void;
  label: string;
}
export default function QuickSelect<
  T extends DefaultMenuSelectOption,
  Multiple extends boolean | undefined = false
>({
  items,
  multiple,
  setSelected,
  selected,
  itemsByName,
  onItemNotFound,
  label,
}: Props<T, Multiple>): JSX.Element {
  const [open, setOpen] = useState(false);
  const [pendingPaste, setPendingPaste] = useState(false);
  const [input, setInput] = useState("");

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

  // `HandleEnter()` handles the enter key press when there is a pending paste in the search bar
  // Since this functionality is currently only used in the gene search bar, we'll be assuming that `itemsByName` is a Map<string, Gene>
  const handleEnter =
    !multiple || !("length" in selected) || onItemNotFound === undefined
      ? noop
      : (event: React.KeyboardEvent<HTMLInputElement>) => {
          if (event.key === "Enter" && pendingPaste) {
            event.preventDefault();
            const newSelected = [...(selected as T[])];
            const parsedPaste = pull(uniq(input.split(/[ ,]+/)), "");
            parsedPaste.map((item) => {
              const newItem = itemsByName.get(item);
              if (!newItem) {
                onItemNotFound(item);
              } else if (!newSelected.includes(newItem))
                newSelected.push(newItem);
            });
            setPendingPaste(false);
            setOpen(false);
            return setSelected(newSelected as Value<T, Multiple>);
          }
        };

  const handlePaste = () => {
    setPendingPaste(true);
  };

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
    newValue: T[] | T | null
  ) => {
    return setSelected(newValue as Value<T, Multiple>);
  };

  const handleClick = () => {
    setOpen(true);
  };

  const ref = useRef(null);

  return (
    <>
      <ButtonWrapper>
        <Label>{label}</Label>
        <StyledIconButton ref={ref} onClick={handleClick} sdsType="primary">
          <Icon sdsIcon="plusCircle" sdsSize="s" sdsType="iconButton" />
        </StyledIconButton>
      </ButtonWrapper>
      <Popper open={open} className={classes.popper} anchorEl={ref.current}>
        <MenuSelect
          open
          search
          onClose={handleClose}
          multiple={"length" in selected}
          classes={{
            paper: classes.paper,
            popperDisablePortal: classes.popperDisablePortal,
          }}
          value={selected}
          onChange={handleChange}
          disableCloseOnSelect
          disableListWrap
          onKeyDownCapture={multiple ? handleEnter : undefined}
          options={items}
          ListboxComponent={
            ListboxComponent as React.ComponentType<
              React.HTMLAttributes<HTMLElement>
            >
          }
          renderOption={(option) => option.name}
          onPaste={handlePaste}
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
    </>
  );
}
