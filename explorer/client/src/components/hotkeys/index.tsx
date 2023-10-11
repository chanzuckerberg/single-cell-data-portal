import { useHotkeys, InputGroup } from "@blueprintjs/core";
import React, { FC, useMemo } from "react";
import { connect } from "react-redux";
import { AppDispatch } from "../../reducers";
import { track } from "../../analytics";
import { EVENTS } from "../../analytics/events";
import * as globals from "../../globals";

interface DispatchProps {
  undo: () => void;
  redo: () => void;
}
type Props = DispatchProps;

const mapDispatchToProps = (dispatch: AppDispatch): DispatchProps => ({
  undo: () => {
    track(EVENTS.EXPLORER_UNDO_BUTTON_CLICKED);
    dispatch({ type: "@@undoable/undo" });
  },
  redo: () => {
    track(EVENTS.EXPLORER_REDO_BUTTON_CLICKED);
    dispatch({ type: "@@undoable/redo" });
  },
});

const GlobalHotkeys: FC<Props> = ({ undo, redo }) => {
  const hotkeys = useMemo(
    () => [
      {
        combo: globals.isMac ? "CMD+Z" : "CTRL+Z",
        global: true,
        label: "Undo.",
        onKeyDown: () => {
          undo();
        },
      },
      {
        combo: globals.isMac ? "CMD+SHIFT+Z" : "CTRL+SHIFT+Z",
        global: true,
        label: "Redo.",
        onKeyDown: () => {
          redo();
        },
      },
    ],
    [],
  );
  const { handleKeyDown, handleKeyUp } = useHotkeys(hotkeys);

  return (
    <div
      role="tab"
      tabIndex={0}
      onKeyDown={handleKeyDown}
      onKeyUp={handleKeyUp}
      style={{
        display: "none",
      }}
    >
      <InputGroup />
    </div>
  );
};

export default connect(null, mapDispatchToProps)(GlobalHotkeys);
