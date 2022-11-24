import { IconNames } from "@blueprintjs/icons";
import { SetSearchValueFn } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import { ClearIconButton, ViewSearch } from "./style";
import { Icon } from "czifui";
import { useRef } from "react";

interface Props {
  className?: string;
  searchValue: string;
  setSearchValue: SetSearchValueFn;
}

export default function FilterSearch({
  className,
  searchValue,
  setSearchValue,
}: Props): JSX.Element {
  const inputRef = useRef<HTMLInputElement>(null);

  /**
   * Clears search value and sets focus on the search input element.
   */
  const onClearSearch = (): void => {
    setSearchValue("");
    inputRef?.current?.focus();
  };

  return (
    <ViewSearch
      autoFocus
      className={className}
      inputRef={inputRef}
      leftIcon={IconNames.SEARCH}
      onChange={(changeEvent) => setSearchValue(changeEvent.target.value)}
      placeholder="Search"
      rightElement={
        searchValue ? (
          <ClearIconButton
            onClick={onClearSearch}
            {...{
              // TODO(cc) move this back to explicit prop={value} after upgrading SDS to enable type checking again
              sdsSize: "small",
            }}
            sdsType="secondary"
          >
            <Icon sdsIcon="xMarkCircle" sdsSize="s" sdsType="iconButton" />
          </ClearIconButton>
        ) : undefined
      }
      value={searchValue}
    />
  );
}
