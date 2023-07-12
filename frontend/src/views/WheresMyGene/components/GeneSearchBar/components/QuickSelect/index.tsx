import { AutocompleteRenderOptionState } from "@mui/material/Autocomplete";
import {
  AutocompleteCloseReason,
  AutocompleteInputChangeReason,
  AutocompleteValue,
} from "@mui/material/useAutocomplete";
import makeStyles from "@mui/styles/makeStyles";
import {
  DefaultMenuSelectOption,
  DropdownPaper,
  DropdownPopper,
  getColors,
  getCorners,
  getShadows,
  Icon,
  MenuSelect,
  SDSTheme,
} from "@czi-sds/components";
import { pull, uniq } from "lodash";

import React, {
  createContext,
  HTMLAttributes,
  MouseEvent,
  ReactChild,
  SyntheticEvent,
  useRef,
  useState,
} from "react";
import { FixedSizeList, ListChildComponentProps } from "react-window";
import { track } from "src/common/analytics";
import { EVENTS } from "src/common/analytics/events";
import { noop } from "src/common/constants/utils";
import { Label } from "../../style";
import {
  ButtonWrapper,
  StyledButtonText,
  StyledMenuItem,
  StyledPopper,
  StyledSelectButton,
} from "./style";

const MAX_ITEMS_TO_SHOW = 9.5;
const LISTBOX_ITEM_HEIGHT_PX = 48;

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

interface ListboxProps {
  children: ReactChild;
}

const ListboxComponent = React.forwardRef<HTMLDivElement, ListboxProps>(
  function ListboxComponent(props: ListboxProps, ref) {
    const { children, ...other } = props;

    const itemData = React.Children.toArray(children);
    const itemCount = itemData.length;

    const height = Math.min(
      LISTBOX_ITEM_HEIGHT_PX * itemCount,
      LISTBOX_ITEM_HEIGHT_PX * MAX_ITEMS_TO_SHOW
    );

    return (
      <div ref={ref}>
        <ListBoxContext.Provider value={other}>
          <FixedSizeList
            height={height}
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

interface Props<T, Multiple> {
  items: T[];
  multiple?: Multiple;
  setSelected: (
    selected: AutocompleteValue<T, Multiple, undefined, undefined>
  ) => void;
  selected: AutocompleteValue<T, Multiple, undefined, undefined>;
  itemsByName: Map<string, T>;
  onItemNotFound?: (item: string) => void;
  label: string;
  text: string;
  dataTestId: string;
  placeholder?: string;
  analyticsEvent: EVENTS;
  isLoading: boolean;
}
export default function QuickSelect<
  T extends DefaultMenuSelectOption,
  Multiple extends boolean | undefined = false
>({
  items,
  multiple,
  setSelected,
  selected,
  /**
   * name is lowercase for case insensitive CSV paste
   */
  itemsByName,
  onItemNotFound,
  label,
  text,
  dataTestId,
  placeholder,
  analyticsEvent,
  isLoading,
}: Props<T, Multiple>): JSX.Element {
  const [open, setOpen] = useState(false);
  const [input, setInput] = useState("");
  const [hasComma, setHasComma] = useState(false);
  const [itemOptions, setItemOptions] = useState<T[]>(items);

  const useStyles = makeStyles((theme: SDSTheme) => {
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
        zIndex: 99, // The x axis wrapper is set at 2
      },
      popperDisablePortal: {
        position: "relative",
        width: "100% !important",
      },
    };
  });

  const classes = useStyles();

  const toggleDropdownAndUpdateItemOptions = () => {
    setOpen((open) => !open);
    const itemOptionsWithSelectedOnTop = [
      ...new Set([...selectedAsArray, ...items]),
    ];
    setItemOptions(itemOptionsWithSelectedOnTop);
  };

  // `HandleEnter()` handles the enter key press when we detect comma(",") in the search bar
  // Since this functionality is currently only used in the gene search bar, we'll be assuming that `itemsByName` is a Map<string, Gene>.
  // NOTE that `itemsByName` the key is lowercase!
  const handleEnter =
    !multiple || !Array.isArray(selected) || onItemNotFound === undefined
      ? noop
      : (event: React.KeyboardEvent<HTMLInputElement>) => {
          if (event.key === "Enter" && hasComma) {
            event.preventDefault();
            const newSelected = [...(selected as T[])];
            const parsedPaste = pull(uniq(input.split(/[ ,]+/)), "");

            parsedPaste.map((item) => {
              const newItem = itemsByName.get(item.toLowerCase());
              if (!newItem) {
                onItemNotFound(item);
              } else if (!newSelected.includes(newItem))
                newSelected.push(newItem);
            });

            toggleDropdownAndUpdateItemOptions();

            return setSelected(
              newSelected as AutocompleteValue<
                T,
                Multiple,
                undefined,
                undefined
              >
            );
          }
        };

  const handleClose = (
    e: SyntheticEvent<Element, Event>,
    reason: AutocompleteCloseReason
  ) => {
    if (reason === "toggleInput") {
      return;
    }
    const { nativeEvent } = e;

    if (
      (nativeEvent instanceof FocusEvent &&
        nativeEvent.relatedTarget instanceof Element &&
        nativeEvent.relatedTarget?.id !== `${dataTestId}-id`) ||
      (nativeEvent instanceof FocusEvent &&
        !(nativeEvent.relatedTarget instanceof Element)) ||
      !(nativeEvent instanceof FocusEvent)
    ) {
      toggleDropdownAndUpdateItemOptions();
    }

    setInput("");
  };
  const handleChange = (
    _: SyntheticEvent<Element, Event>,
    newValue: DefaultMenuSelectOption | DefaultMenuSelectOption[] | null
  ) => {
    return setSelected(
      newValue as AutocompleteValue<T, Multiple, undefined, undefined>
    );
  };

  const handleClick = () => {
    toggleDropdownAndUpdateItemOptions();
  };

  const ref = useRef(null);

  const handleInputChange = (
    _: SyntheticEvent<Element, Event>,
    value: string,
    reason: AutocompleteInputChangeReason
  ) => {
    if (reason === "reset") {
      return;
    }

    setInput(value);

    if (value.includes(",")) {
      if (hasComma) return;
      setHasComma(true);
    } else {
      if (!hasComma) return;
      setHasComma(false);
    }
  };
  const selectedAsArray = Array.isArray(selected) ? selected : [selected];

  return (
    <>
      <ButtonWrapper>
        <Label>{label}</Label>
        <StyledSelectButton
          disabled={isLoading}
          id={`${dataTestId}-id`}
          data-testid={dataTestId}
          ref={ref}
          onClick={handleClick}
          sdsType="primary"
          sdsSize="medium"
          sdsStyle="square"
        >
          <StyledButtonText>
            <Icon sdsIcon="plus" sdsSize="s" sdsType="button" />
          </StyledButtonText>
          {text}
        </StyledSelectButton>
      </ButtonWrapper>
      <StyledPopper
        open={open}
        className={classes.popper}
        anchorEl={ref.current}
      >
        <MenuSelect
          open
          PopperComponent={DropdownPopper}
          PaperComponent={DropdownPaper}
          search
          onClose={handleClose}
          multiple={multiple}
          classes={{
            paper: classes.paper,
            popperDisablePortal: classes.popperDisablePortal,
          }}
          value={selected}
          onChange={handleChange}
          disableCloseOnSelect
          disableListWrap
          onKeyDownCapture={multiple ? handleEnter : undefined}
          options={itemOptions}
          ListboxComponent={
            ListboxComponent as React.ComponentType<
              React.HTMLAttributes<HTMLElement>
            >
          }
          renderOption={renderOption}
          InputBaseProps={{
            placeholder,
          }}
          inputValue={input}
          onInputChange={handleInputChange}
          noOptionsText={
            hasComma
              ? "You can add multiple genes using a comma-separated list. Press enter to add."
              : "No options"
          }
        />
      </StyledPopper>
    </>
  );

  function renderOption(
    /**
     * (thuang): MUI passes their own event handlers via `optionProps`, so if we
     * need our own event handlers, we need to make sure to call their handlers
     * too
     */
    optionProps: HTMLAttributes<HTMLLIElement>,
    option: T,
    { selected }: AutocompleteRenderOptionState
  ) {
    return (
      <StyledMenuItem
        {...{ component: "li" }}
        {...optionProps}
        isMultiSelect={multiple}
        selected={selected}
        onClick={onClick}
      >
        {option.name}
      </StyledMenuItem>
    );

    function onClick(event: MouseEvent<HTMLLIElement>): void {
      optionProps.onClick && optionProps.onClick(event);

      // (thuang): Only track select, not deselect
      if (selected) return;

      track(analyticsEvent, { payload: option.name });
    }
  }
}
