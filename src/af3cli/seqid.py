from string import ascii_uppercase


def num_to_letters(num: int) -> str:
    """
    Convert a positive integer to its corresponding uppercase
    alphabetical representation.

    Parameters
    ----------
    num : int
        A positive integer to convert to an alphabetical
        representation. Must be greater than 0.

    Returns
    -------
    str
        The corresponding uppercase alphabetical representation of
        the given integer.

    Raises
    ------
    ValueError
        If the input `num` is less than or equal to 0, an exception
        is raised indicating an invalid input.

    Notes
    -----
    The mapping is consecutive, with 1 corresponding to 'A', 2 to 'B', ...,
    and 27 to 'AA'.
    """
    if num <= 0:
        raise ValueError("Sequence ID count must be greater than 0")
    result = ""
    while num > 0:
        num -= 1
        result = ascii_uppercase[num % 26] + result
        num //= 26
    return result


class IDRecord(object):
    """
    Manages a sequence ID record to automatically handle sequence IDs
    in the AlphaFold3 input file.

    Attributes
    ----------
    _seq_id : list of str or None
        The stored sequence IDs.
    _num : int or None
        The number of ligand sequences, default is 1. This value
        will be overwritten if `seq_id` is larger.
    """
    def __init__(self, num: int = 1, seq_id: list[str] | None = None):
        self._seq_id: list[str] | None = seq_id
        self._tmp_seq_id: list[str] = []
        self._num: int = num
        self._sanitize_seq_id(seq_id)
        self._sanitize_num(num)

    def _sanitize_seq_id(self, seq_id: list[str] | str | None) -> None:
        """
        Sanitizes the sequence identifier.

        This method processes the provided sequence identifier `seq_id`. If `seq_id`
        is None or contains no elements, the internal `_seq_id` attribute is set
        to None. Otherwise, `_seq_id` is updated with the given `seq_id`.

        Parameters
        ----------
        seq_id : list of str or None
            A list of sequence identifiers. If None or an empty list is
            provided, the internal `_seq_id` attribute is set to None; otherwise,
            it is assigned the provided list.
        """
        if isinstance(seq_id, str) and len(seq_id) > 0:
            self._seq_id = [seq_id]
            return
        if seq_id is None or len(seq_id) == 0:
            self._seq_id = None
        else:
            self._seq_id = seq_id

    def _sanitize_num(self, num: int = 1) -> None:
        """
        Sanitizes and updates the `_num` attribute based on given input or the
        length of `_seq_id`.

        Parameters
        ----------
        num : int, optional
            An integer value indicating the desired number. Defaults to 1.
        """
        if num < 1:
            raise ValueError("Sequence ID count must be greater than 0")
        if self._seq_id is None or len(self._seq_id) <= num:
            self._num = num
        else:
            self._num = len(self._seq_id)

    @property
    def num(self) -> int:
        return self._num

    @num.setter
    def num(self, num: int) -> None:
        self._sanitize_num(num)

    def set_temporary_id(self, seq_id: list[str]) -> None:
        """
        Sets a temporary sequence identifier for the object.

        Parameters
        ----------
        seq_id : list of str
            The list of string identifiers to assign as the temporary
            sequence identifier.
        """
        self._tmp_seq_id = seq_id

    def get_temporary_id(self) -> list[str]:
        return self._tmp_seq_id

    def required_tmp_id_count(self):
        if self._seq_id is None:
            return self._num
        return self._num - len(self._seq_id)

    def get_id(self) -> list[str] | None:
        return self._seq_id

    def set_id(self, seq_id: list[str] | str | None) -> None:
        """
        Set a new sequence identifier for the object. If an empty list is
        passed as the sequence identifier, the identifier will be set to None.

        Parameters
        ----------
        seq_id : list of str or None
            A list of string sequence identifiers to be associated with the
            object. If an empty list is passed, the sequence identifier is set
            to None.
        """
        self._sanitize_seq_id(seq_id)
        self._sanitize_num(self._num)

    def get_full_id_list(self) -> list[str]:
        """
        Gets the complete list of IDs by combining the explicitly set IDs
        and temporary IDs.

        Returns
        -------
        list of str
            A combined list of IDs from the primary and temporary list, or just
            the temporary list if the list of explicit IDs is None.
        """
        if self._seq_id is None:
            return self._tmp_seq_id
        return self._seq_id + self._tmp_seq_id

    def remove_id(self) -> None:
        self.set_id(None)

    def clear_temporary_id(self) -> None:
        self._tmp_seq_id = []


class IDRegister(object):
    """
    Manages the registration and generation of unique sequence IDs.

    This class ensures that sequence IDs are unique and allows for the
    registration of existing IDs as well as the generation of new, unique
    IDs based on an internal counter. It is primarily designed to prevent
    duplication and collisions of IDs in the AlphaFold3 input file.

    Attributes
    ----------
    _count : int
        The internal counter used for generating unique sequence IDs.
    _registered_ids : set
        A set that stores all the registered sequence IDs to ensure
        uniqueness.
    """
    def __init__(self):
        self._count = 0
        self._registered_ids = set()

    def register(self, seq_id: str) -> None:
        """
        Registers a new sequence ID, ensuring it is unique.

        Parameters
        ----------
        seq_id : str
            The new sequence ID to register.

        Raises
        ------
        ValueError
            If the sequence ID has already been registered.

        """
        if seq_id in self._registered_ids:
            raise ValueError(
                f"Sequence ID {seq_id} has already been registered"
            )
        self._registered_ids.add(seq_id)

    def generate(self) -> str:
        """
        Generates a new unique sequence ID.

        Returns
        -------
        str
            A new unique sequence ID.
        """
        self._count += 1
        while True:
            seq_id = num_to_letters(self._count)
            if seq_id not in self._registered_ids:
                return seq_id
            self._count += 1

    def reset(self) -> None:
        """
        Resets the `IDRegister` object.

        Notes
        -----
        The ids and register states of corresponding sequence objects are
        not affected.
        """
        self._count = 0
        self._registered_ids = set()
